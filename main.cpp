#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using namespace std;

vector<string> readFile(string fileName) {
    vector<string> lines;
    fstream myfile (fileName);
    string line;
    while (getline(myfile, line))
        lines.push_back(line);
    myfile.close();
    lines.erase(lines.begin());
    return lines;
}

vector<string> split(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<string> elems;
    while (getline(ss, item, delim))
        elems.push_back(item);
    return elems;
}

class Branch {
    int nodeI, nodeJ, branchK, type; // type = 1 for R, type = 2 for E, type = 3 for I
    double value;
public:
    Branch(int nodeI, int nodeJ, int branchK, int type, double value) : nodeI(nodeI), nodeJ(nodeJ), branchK(branchK),
                                                                        type(type), value(value) {}
    int getNodeI() const { return nodeI; }
    int getNodeJ() const { return nodeJ; }
    int getBranchK() const { return branchK; }
    int getType() const { return type; }
    double getValue() const { return value; }

    void setNodeI(int nodeI) { Branch::nodeI = nodeI; }
    void setNodeJ(int nodeJ) { Branch::nodeJ = nodeJ; }
    void setBranchK(int branchK) { Branch::branchK = branchK; }
    void setType(int type) { Branch::type = type; }
    void setValue(double value) { Branch::value = value; }
};

class Circuit {
    int noOfNodes, refNode, n, m;
    bool solved;
    vector<double> nodesVoltages, volSourcesCurrents;
    vector<Branch> branches, noRefBranches;

    int noOfVolSources() {
        int n = 0;
        for (Branch b:branches) {
            if (b.getType() == 2) n++;
        }
        return n;
    }
    MatrixXd admittances() {
        MatrixXd G = MatrixXd::Zero(n, n);
        for (Branch b:branches) {
            if (b.getType() != 1) continue;
            int i = b.getNodeI(), j = b.getNodeJ();
            j = (j < refNode) ? j : j - 1;
            i = (i < refNode) ? i : i - 1;
            double v = b.getValue();
            if (i != 0) {
                G(i-1, j-1) -= 1/v;
                G(j-1, i-1) -= 1/v;
                G(i-1, i-1) += 1/v;
            }
            G(j-1, j-1) += 1/v;
        }
        return G;
    }
    MatrixXd volSourcesConnections() {
        MatrixXd B = MatrixXd::Zero(n, m);
        int k = 0;
        for (Branch b:branches) {
            if (b.getType() != 2) continue;
            int i = b.getNodeI(), j = b.getNodeJ(), v = -1;
            j = (j < refNode) ? j : j - 1;
            i = (i < refNode) ? i : i - 1;
            if(b.getValue() > 0) v = 1;
            if (i != 0)
                B(i-1, k) = -v;
            B(j-1, k) = v;
            k++;
        }
        return B;
    }
    MatrixXd currentSources() {
        MatrixXd I = MatrixXd::Zero(n, 1);
        for (Branch b:branches) {
            if (b.getType() != 3) continue;
            int i = b.getNodeI(), j = b.getNodeJ();
            j = (j < refNode) ? j : j - 1;
            i = (i < refNode) ? i : i - 1;
            double v = b.getValue();
            if (i != 0)
                I(i-1, 0) -= v;
            I(j-1, 0) += v;
        }
        return I;
    }
    MatrixXd voltageSources() {
        MatrixXd e = MatrixXd::Zero(m, 1);
        int k = 0;
        for (Branch b:branches) {
            if (b.getType() == 2) {
                e(k, 0) = abs(b.getValue());
                k++;
            }
        }
        return e;
    }

    void countK(vector<Branch> &vecB) {
        for (int x(0); x < vecB.size(); x++) {
            int nodeIx = vecB[x].getNodeI(), nodeJx = vecB[x].getNodeJ();
            int k = 1;
            bool a = false;
            for (int y(x + 1); y < vecB.size(); y++) {
                int nodeIy = vecB[y].getNodeI(), nodeJy = vecB[y].getNodeJ();
                if(nodeIx == nodeIy && nodeJx == nodeJy && vecB[y].getBranchK() == 0) {
                    vecB[y].setBranchK(++k);
                    a = true;
                }
            }
            if(a) vecB[x].setBranchK(1);
        }
    }

public:
    Circuit(int noOfNodes) : noOfNodes(noOfNodes) { solved = false; }

    Circuit(string fileName) {
        vector<string> lines = readFile(fileName), types, values, is, js;
        vector<Branch> vecB;
        int num = 1;

        // filling vectors "types", "values", "is" and "js"
        for(string s:lines) {
            vector<string> parts = split(s, ' ');
            if(s == "" || (parts[0] != "r" && parts[0] != "v" && parts[0] != "i" && parts[0] != "w")) continue;
            types.push_back(parts[0]);
            is.push_back(parts[1] + "," + parts[2]);
            js.push_back(parts[3] + "," + parts[4]);
            string value = (parts[0] == "v") ? parts[8] : parts[parts.size() - 1];
            values.push_back(value);
        }
        // making one point for same potential
        for (int k(0); k < lines.size(); k++) {
            if(types[k] == "w") {
                for (int l(0); l < lines.size(); l++) {
                    string wireI = is[k];
                    if(is[l] == wireI) is[l] = js[k];
                    if(js[l] == wireI) js[l] = js[k];
                }
            }
        }
        // erasing wires
        for (int k(0); k < types.size(); k++) {
            if(types[k] == "w") {
                types.erase(types.begin() + k);
                is.erase(is.begin() + k);
                js.erase(js.begin() + k);
                values.erase(values.begin() + k);
                k--;
            }
        }
        // replacing positions with numbers beginning from 1
        for (int k(0); k < types.size(); k++) {
            if(types[k] != "w") {
                string compI = is[k], compJ = js[k];
                if(compI.find(',') != string::npos) {
                    for (int l(0); l < types.size(); l++) {
                        if(is[l] == compI) is[l] = to_string(num);
                        if(js[l] == compI) js[l] = to_string(num);
                    }
                    num++;
                }
                if(compJ.find(',') != string::npos) {
                    for (int l(0); l < types.size(); l++) {
                        if(is[l] == compJ) is[l] = to_string(num);
                        if(js[l] == compJ) js[l] = to_string(num);
                    }
                    num++;
                }
            }
        }
        for (int k(0); k < types.size(); k++) {
            int type, i, j;
            double value;
            type = (types[k] == "r") ? 1 : ((types[k] == "v") ? 2 : 3);
            i = stoi(is[k]);
            j = stoi(js[k]);
            value = stod(values[k]);

            if(i > j) {
                int temp = i;
                i = j;
                j = temp;
                if(type != 1) value = -value;
            }
            Branch branch(i, j, 0, type, value);
            vecB.push_back(branch);
        }
        noOfNodes = num - 1;
        countK(vecB);
        setBranches(vecB);
        solved = false;
    }

    int getNoOfNodes() const { return noOfNodes; }
    int getRefNode() const { return refNode; }
    const vector<Branch> &getBranches() const { return branches; }
    const vector<double> &getNodesVoltages() const { return nodesVoltages; }
    const vector<double> &getVolSourcesCurrents() const { return volSourcesCurrents; }

    void setNoOfNodes(int noOfNodes) { Circuit::noOfNodes = noOfNodes; solved = false; }
    void setRefNode(int refNode) {
        if(refNode < 1 || refNode > noOfNodes) throw "Reference node out of boundaries";
        Circuit::refNode = refNode;
        for (int i(0); i < branches.size(); i++) {
            if (branches[i].getNodeI() == refNode) {
                branches[i].setNodeI(0);
            } else if (branches[i].getNodeJ() == refNode) {
                branches[i].setNodeJ(branches[i].getNodeI());
                branches[i].setNodeI(0);
                if(branches[i].getType() != 1) branches[i].setValue(-branches[i].getValue());
            }
        }
        solved = false;
    }
    void setBranches(const vector<Branch> &branches) {
        Circuit::branches = branches;
        noRefBranches = branches;
        solved = false;
    }

    void solve() {
        n = noOfNodes - 1, m = noOfVolSources();
        MatrixXd G(n, n), B(n, m), D = MatrixXd::Zero(m, m), A(n+m, n+m), i(n, 1), e(m, 1), b(n+m, 1), x;
        G = admittances();
        B = volSourcesConnections();
        A << G,             B,
             B.transpose(), D;
        i = currentSources();
        e = voltageSources();
        b << i, e;

        x = A.inverse() * b;
        MatrixXd vn(noOfNodes, 1), iv = x.block(n, 0, m, 1);
        vn << x.block(0, 0, refNode-1, 1), 0, x.block(refNode-1, 0, noOfNodes-refNode, 1); // inserting 0 at refNode index
        vector<double> vecVn(vn.data(), vn.data() + vn.size());
        vector<double> vecIv(iv.data(), iv.data() + iv.size());
        nodesVoltages = vecVn;
        volSourcesCurrents = vecIv;
        solved = true;
    }
    vector<double> getBranchesCurrents() {
        if(!solved) throw "Cicuit is not solved yet";
        vector<double> currents;
        int vx = 0;
        for (Branch b:noRefBranches) {
            if (b.getType() == 3)
                currents.push_back(b.getValue());
            else if (b.getType() == 2) {
                if(b.getValue() > 0) currents.push_back(-volSourcesCurrents[vx++]);
                else                 currents.push_back(volSourcesCurrents[vx++]);
            }
            else if (b.getType() == 1) {
                double vi = nodesVoltages[b.getNodeI() - 1], vj = nodesVoltages[b.getNodeJ() - 1];
                currents.push_back((vi - vj) / b.getValue());
            }
        }
        return currents;
    }
    void printBranchesCurrents() {
        if(!solved) throw "Circuit is not solved yet";
        vector<double> currents = getBranchesCurrents();
        for (int k(0); k < currents.size(); k++) {
            cout << "I_" << noRefBranches[k].getNodeI() << "_" << noRefBranches[k].getNodeJ();
            if (noRefBranches[k].getBranchK() != 0) cout << "_" << noRefBranches[k].getBranchK();
            printf(" = %.3lfmA \n", currents[k] * 1000);
        }
    }
};

void program() {
    cout << "\n--------------------------------------------------------------\n"
            "DCCalculator - Aplikacija za proračun DC mreže korištenjem MNČ\n"
            "--------------------------------------------------------------\n\n"
            "Prvo se unese broj čvorova (dijelovi  mreže  koji  nemaju isti\n"
            "potencijal), onda  se unese broj grana između svaka dva čvora,\n"
            "nakon toga  se unese šta svaka grana sadrži i vrijednost toga,\n"
            "te unese se broj referentnog čvora.\n"
            "Kliknite enter za nastavak...";
    cin.ignore();
    int n, m, k, kTemp, type, refNode;
    double value;
    vector<Branch> branches;
    cout << "\nBroj čvorova u mreži: ";
    cin >> n;
    for(int i(1); i <= n - 1; i++) {
        for(int j(i + 1); j <= n; j++) {
            cout << "\nBroj grana između " << i << ". i " << j << ". čvor (0 ako nisu povezani): ";
            cin >> m;
            for(int k(1); k <= m; k++) {
                if(m == 1)  cout << "Grana između " << i << ". i " << j << ". čvor sadrži (R - 1, E - 2, Is - 3): ";
                else cout << k << ". grana između " << i << ". i " << j << ". čvor sadrži (R - 1, E - 2, Is - 3): ";
                cin >> type;
                if (type == 2) cout << "Vrijednosti naponskog izvora (pomnožiti sa -1 ako je + na strani čvora " << i <<"): ";
                else if (type == 3) cout << "Vrijednosti strujnog izvora (pomnožiti sa -1 ako strelica ide ka čvoru " << i <<"): ";
                else cout << "Vrijednosti otpornika: ";
                cin >> value;
                if(m == 1) kTemp = 0;
                else kTemp = k;
                Branch branch(i,j,kTemp,type,value);
                branches.push_back(branch);
            }
        }
    }
    cout << "\nReferentni čvor: ";
    cin >> refNode;
    cout<<endl;

    Circuit cir(n);
    cir.setBranches(branches);
    cir.setRefNode(refNode);
    cir.solve();
    cir.printBranchesCurrents();
}

int main() {
    program();
}
