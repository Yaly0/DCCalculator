#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using namespace std;

vector<string> readFile(const string &fileName) {
    vector<string> lines;
    fstream file(fileName);
    string line;
    while (getline(file, line))
        lines.push_back(line);
    file.close();
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
    int nodeI, nodeJ, branchK, type;
    // type = 1 for R:resistor, 2 - E:volSource, 3 - Is:currSource, 4 - Uv:voltmeter, 5 - Ia:ammeter,
    // 6 - calculated ammeter, 7 - Wv: voltage of wattmeter, 8 - Wa: current of wattmeter, 9 - calculated Wa
    double value;
public:
    Branch(int nodeI, int nodeJ, int branchK, int type, double value) : nodeI(nodeI), nodeJ(nodeJ), branchK(branchK),
                                                                        type(type), value(value) {}
    int getNodeI() const { return nodeI; }
    int getNodeJ() const { return nodeJ; }
    int getBranchK() const { return branchK; }
    int getType() const { return type; }
    double getValue() const { return value; }

    void setNodeI(int nI) { Branch::nodeI = nI; }
    void setNodeJ(int nJ) { Branch::nodeJ = nJ; }
    void setBranchK(int bK) { Branch::branchK = bK; }
    void setType(int t) { Branch::type = t; }
    void setValue(double v) { Branch::value = v; }
};

ostream &operator<<(ostream &Str, Branch const &b) {
    Str << to_string(b.getType()) << " " << to_string(b.getNodeI()) << ","
        << to_string(b.getNodeJ()) << " "<< to_string(b.getValue());
    return Str;
}

void appendToFile(const string& fileName, Branch b) {
    ofstream file(fileName, ios::app);
    int t = b.getType(), i = b.getNodeI(), j = b.getNodeJ();
    string type;
    if( t == 1) type = "r"; else if( t == 2) type = "v"; else if( t == 3) type = "i"; else if( t == 4) type = "p";
    else if( t == 5) type = "370"; else if( t == 7) type = "p420"; else if( t == 8) type = "i420";
    if (type == "v")
        file << type << " " << i << " " << i << " " << j << " " << j << " 0 0 0 " << b.getValue() << "\n";
    else
        file << type << " " << i << " " << i << " " << j << " " << j << " " << b.getValue() << "\n";
    file.close();
}

class Circuit {
    int noOfNodes{}, refNode{}, n{}, m{};
    bool solved;
    vector<double> nodesVoltages, volSourcesCurrents;
    vector<Branch> branches, noRefBranches;

    int noOfVolSources() {
        int i = 0;
        for (Branch b: branches) {
            if (b.getType() == 2) i++;
        }
        return i;
    }
    MatrixXd admittances() {
        MatrixXd G = MatrixXd::Zero(n, n);
        for (Branch b: branches) {
            if (b.getType() != 1) continue;
            int i = b.getNodeI(), j = b.getNodeJ();
            j = (j < refNode) ? j : j - 1;
            i = (i < refNode) ? i : i - 1;
            double v = b.getValue();
            if (i == 0) {
                G(j - 1, j - 1) += 1 / v;
                continue;
            }
            if (j == 0) {
                G(i - 1, i - 1) += 1 / v;
                continue;
            }
            G(i - 1, j - 1) -= 1 / v;
            G(j - 1, i - 1) -= 1 / v;
            G(i - 1, i - 1) += 1 / v;
            G(j - 1, j - 1) += 1 / v;
        }
        return G;
    }
    MatrixXd volSourcesConnections() {
        MatrixXd B = MatrixXd::Zero(n, m);
        int k = 0;
        for (Branch b: branches) {
            if (b.getType() != 2) continue;
            int i = b.getNodeI(), j = b.getNodeJ(), v = -1;
            j = (j < refNode) ? j : j - 1;
            i = (i < refNode) ? i : i - 1;
            if (b.getValue() > 0) v = 1;
            if (i != 0) B(i - 1, k) = -v;
            if (j != 0) B(j - 1, k) = v;
            k++;
        }
        return B;
    }
    MatrixXd currentSources() {
        MatrixXd I = MatrixXd::Zero(n, 1);
        for (Branch b: branches) {
            if (b.getType() != 3) continue;
            int i = b.getNodeI(), j = b.getNodeJ();
            j = (j < refNode) ? j : j - 1;
            i = (i < refNode) ? i : i - 1;
            double v = b.getValue();
            if (i != 0) I(i - 1, 0) -= v;
            if (j != 0) I(j - 1, 0) += v;
        }
        return I;
    }
    MatrixXd voltageSources() {
        MatrixXd e = MatrixXd::Zero(m, 1);
        int k = 0;
        for (Branch b: branches) {
            if (b.getType() == 2) {
                e(k, 0) = abs(b.getValue());
                k++;
            }
        }
        return e;
    }

    static void countK(vector<Branch> &vecB) {
        for (int x(0); x < vecB.size(); x++) {
            int nodeIx = vecB[x].getNodeI(), nodeJx = vecB[x].getNodeJ();
            int k = 1;
            bool a = false;
            for (int y(x + 1); y < vecB.size(); y++) {
                int nodeIy = vecB[y].getNodeI(), nodeJy = vecB[y].getNodeJ();
                if ((nodeIx == nodeIy && nodeJx == nodeJy || nodeIx == nodeJy && nodeJx == nodeIy) &&
                    vecB[y].getBranchK() == 0) {
                    vecB[y].setBranchK(++k);
                    a = true;
                }
            }
            if (a) vecB[x].setBranchK(1);
        }
    }
    vector<Branch> branchesFromTxt(const string& fileName, bool withAmmeters) {
        vector<string> lines = readFile(fileName), types, values, is, js;
        vector<Branch> vecB;
        int num = 1;

        // filling vectors "types", "values", "is" and "js"
        for (int k(0); k < lines.size(); k++) {
            vector<string> parts = split(lines[k], ' ');
            if (lines[k].empty() || (parts[0] != "r" && parts[0] != "v" && parts[0] != "i" &&
                                     parts[0] != "w" && parts[0] != "p" && parts[0] != "370" &&
                                     parts[0] != "420" && parts[0] != "p420" && parts[0] != "i420"))
                continue;
            if (parts[0] == "420") { // dividing wattmeter into ammeter and voltmeter and a wire
                int vX, vY, wX, wY;
                if (parts[2] == parts[4]) { // horizontal
                    vX = stoi(parts[1]);
                    vY = stoi(parts[2]) + stoi(parts[6]);
                    wX = stoi(parts[3]);
                    wY = vY;
                } else if (parts[1] == parts[3]) {  // vertical
                    vX = stoi(parts[1]) + stoi(parts[6]);
                    vY = stoi(parts[2]);
                    wX = vX;
                    wY = stoi(parts[4]);
                } else throw logic_error("Wattmeter should be placed vertically or horizontally");
                string vMeter =
                        "p420 " + parts[1] + " " + parts[2] + " " + to_string(vX) + " " + to_string(vY) + " 1 0 0";
                string wire =
                        "w " + to_string(vX) + " " + to_string(vY) + " " + to_string(wX) + " " + to_string(wY) + " 0";
                lines.insert(lines.begin() + k + 1, wire);
                lines.insert(lines.begin() + k + 1, vMeter);
                types.emplace_back("i420");
                goto skip;
            }
            types.push_back(parts[0]);
            skip:
            is.push_back(parts[1] + "," + parts[2]);
            js.push_back(parts[3] + "," + parts[4]);
            string value = (parts[0] == "v") ? parts[8] : parts[parts.size() - 1];
            values.push_back(value);
        }

        if (types.empty()) throw logic_error("Prazan krug");

        string ammeter = "370", wattmeter = "i420";
        if (withAmmeters) {
            ammeter = "skip";
            wattmeter = "skip";
        }

        // making one point for same potential
        for (int k(0); k < types.size(); k++) {
            if (types[k] == "w" || types[k] == ammeter || types[k] == wattmeter) {
                string wireI = is[k], wireJ = js[k];
                for (int l(0); l < types.size(); l++) {
                    if (is[l] == wireI) is[l] = wireJ;
                    if (js[l] == wireI) js[l] = wireJ;
                }
            }
        }
        if (withAmmeters) {
            ammeter = "p";
            wattmeter = "p420";
        }
        // erasing wires
        for (int k(0); k < types.size(); k++) {
            if (types[k] == "w" || types[k] == ammeter || types[k] == wattmeter) {
                types.erase(types.begin() + k);
                is.erase(is.begin() + k);
                js.erase(js.begin() + k);
                values.erase(values.begin() + k);
                k--;
            }
        }
        // replacing positions with numbers beginning from 1
        for (int k(0); k < types.size(); k++) {
            string compI = is[k], compJ = js[k];
            if (compI.find(',') != string::npos) {
                for (int l(0); l < types.size(); l++) {
                    if (is[l] == compI) is[l] = to_string(num);
                    if (js[l] == compI) js[l] = to_string(num);
                }
                num++;
            }
            if (compJ.find(',') != string::npos) {
                for (int l(0); l < types.size(); l++) {
                    if (is[l] == compJ) is[l] = to_string(num);
                    if (js[l] == compJ) js[l] = to_string(num);
                }
                num++;
            }
        }
        for (int k(0); k < types.size(); k++) {
            int type, i, j;
            double value;
            if (types[k] == "r") type = 1; // Resistor
            else if (types[k] == "v") type = 2; // Voltage source
            else if (types[k] == "i") type = 3; // Current source
            else if (types[k] == "p") type = 4; // Voltmeter
            else if (types[k] == "370") type = 5; // Ammeter
            else if (types[k] == "p420") type = 7; // Wattmeter voltage
            else if (types[k] == "i420") type = 8; // Wattmeter current

            i = stoi(is[k]);
            j = stoi(js[k]);
            value = stod(values[k]);
            Branch branch(i, j, 0, type, value);
            vecB.push_back(branch);
        }
        noOfNodes = num - 1;
        countK(vecB);
        return vecB;
    }

public:
    explicit Circuit(int noOfNodes) : noOfNodes(noOfNodes) { solved = false; }
    explicit Circuit(const string &fileName) {
        vector<Branch> vecB = branchesFromTxt((fileName), false);
        setBranches(vecB);
        solved = false;
    }

    void setRefNode(int rNode) {
        if (rNode < 1 || rNode > noOfNodes) throw out_of_range("Reference node out of boundaries");
        Circuit::refNode = rNode;
        for (Branch branch: branches) {
            if (branch.getNodeI() == rNode) {
                branch.setNodeI(0);
            } else if (branch.getNodeJ() == rNode) {
                branch.setNodeJ(branch.getNodeI());
                branch.setNodeI(0);
                if (branch.getType() != 1) branch.setValue(-branch.getValue());
            }
        }
        solved = false;
    }
    void setBranches(const vector<Branch> &vecB) {
        Circuit::branches = vecB;
        noRefBranches = vecB;
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
        if (!solved) throw logic_error("Circuit is not solved yet");
        vector<double> currents;
        int vx = 0;
        for (Branch b: noRefBranches) {
            if (b.getType() == 3)
                currents.push_back(b.getValue());
            else if (b.getType() == 2) {
                if (b.getValue() > 0) currents.push_back(-volSourcesCurrents[vx++]);
                else currents.push_back(volSourcesCurrents[vx++]);
            } else if (b.getType() == 1) {
                double vi = nodesVoltages[b.getNodeI() - 1], vj = nodesVoltages[b.getNodeJ() - 1];
                currents.push_back((vi - vj) / b.getValue());
            }
        }
        return currents;
    }
    void printVoltmeters() {
        if (!solved) throw logic_error("Circuit is not solved yet");
        vector<double> voltmeters;
        for (Branch b: noRefBranches) {
            if (b.getType() == 4) {
                int vi = b.getNodeI(), vj = b.getNodeJ();
                double voltage = nodesVoltages[vi - 1] - nodesVoltages[vj - 1];
                voltmeters.push_back(voltage);
            }
        }
        if (!voltmeters.empty()) cout << "\nVoltmetri:\n";
        for (int k(0); k < voltmeters.size(); k++) {
            cout << "Uv_" << k + 1;
            printf(" = %.3lfV\n", voltmeters[k]);
        }
    }
    void printCurrents() {
        if (!solved) throw logic_error("Circuit is not solved yet");
        vector<Branch> newBranches = branchesFromTxt("falstad.txt", true);
        vector<double> currents = getBranchesCurrents();
        bool hasAmmeters = false;
        for (int k(0); k < newBranches.size(); k++) {
            if (newBranches[k].getType() == 5 || newBranches[k].getType() == 8)
                currents.insert(currents.begin() + k, 0);
        }
        cout << endl;
        int connections, index, sign, ki, kj, li, lj;
        for (int k(0); k < newBranches.size(); k++) { // Ammeters in series with other components
            ki = newBranches[k].getNodeI();
            kj = newBranches[k].getNodeJ();
            if (newBranches[k].getType() == 5 || newBranches[k].getType() == 8) {
                if (newBranches[k].getType() == 5) hasAmmeters = true;
                sign = 1;
                connections = 0;
                for (int l(0); l < newBranches.size(); l++) {
                    li = newBranches[l].getNodeI();
                    lj = newBranches[l].getNodeJ();
                    if (k != l && (li == ki || lj == ki)) {
                        if (li == ki) sign = -1;
                        connections++;
                        index = l;
                    }
                }
                if (connections == 1) {
                    currents[k] = sign * currents[index];
                    newBranches[k].setType(newBranches[k].getType() + 1);
                    continue;
                }
                sign = 1;
                connections = 0;
                for (int l(0); l < newBranches.size(); l++) {
                    li = newBranches[l].getNodeI();
                    lj = newBranches[l].getNodeJ();
                    if (k != l && (li == kj || lj == kj)) {
                        if (lj == kj) sign = -1;
                        connections++;
                        index = l;
                    }
                }
                if (connections == 1) {
                    currents[k] = sign * currents[index];
                    newBranches[k].setType(newBranches[k].getType() + 1);
                    continue;
                }
            }
        }
        double current;
        for (int k(0); k < newBranches.size(); k++) { // Ammeters alone in a wire
            ki = newBranches[k].getNodeI();
            kj = newBranches[k].getNodeJ();
            if (newBranches[k].getType() == 5 || newBranches[k].getType() == 8) {
                current = 0;
                for (int l(0); l < newBranches.size(); l++) {
                    li = newBranches[l].getNodeI();
                    lj = newBranches[l].getNodeJ();
                    if (k != l && (li == ki || lj == ki)) {
                        if (newBranches[l].getType() == 5 || newBranches[l].getType() == 8) goto nodeJ;
                        if (li == ki) current -= currents[l];
                        else current += currents[l];
                    }
                }
                currents[k] = current;
                newBranches[k].setType(newBranches[k].getType() + 1);
                continue;
                nodeJ:
                current = 0;
                for (int l(0); l < newBranches.size(); l++) {
                    li = newBranches[l].getNodeI();
                    lj = newBranches[l].getNodeJ();
                    if (k != l && (li == kj || lj == kj)) {
                        if (newBranches[l].getType() == 5 || newBranches[l].getType() == 8) goto next;
                        if (lj == kj) current += currents[l];
                        else current -= currents[l];
                    }
                }
                currents[k] = current;
                newBranches[k].setType(newBranches[k].getType() + 1);
                next:;
            }
        }

        if (!currents.empty()) cout << "Struje kroz grane:\n";
        for (int k(0); k < currents.size(); k++) { // printing branches' currents
            if (newBranches[k].getType() > 3) continue;
            cout << "I_" << newBranches[k].getNodeI() << "_" << newBranches[k].getNodeJ();
            if (newBranches[k].getBranchK() != 0) cout << "_" << newBranches[k].getBranchK();
            printf(" = %.3lfmA\n", abs(currents[k] * 1000) < 0.0005 ? 0.000 : currents[k] * 1000); // :? operator used to avoid -0
        }

        if (hasAmmeters) cout << "\nAmpermetri:\n";
        int i = 1;
        for (int k(0); k < currents.size(); k++) { // printing ammeters' currents
            if (newBranches[k].getType() == 6) {
                cout << "Ia_" << i++;
                printf(" = %.3lfmA\n", currents[k] * 1000);
            }
        }
        vector<double> wattmetersVoltages;
        for (Branch b: noRefBranches) {
            if (b.getType() == 7) {
                int vi = b.getNodeI(), vj = b.getNodeJ();
                double voltage = nodesVoltages[vi - 1] - nodesVoltages[vj - 1];
                wattmetersVoltages.push_back(voltage);
            }
        }
        if (!wattmetersVoltages.empty()) cout << "\nVatmetri:\n";
        i = 1;
        for (int k(0); k < currents.size(); k++) { // printing wattmeters' powers
            if (newBranches[k].getType() == 9) {
                cout << "Pw_" << i++;
                printf(" = %.3lfmW\n", currents[k] * wattmetersVoltages[i - 2] * 1000);
            }
        }
    }
    void printSolution() {
        if (!solved) throw logic_error("Circuit is not solved yet");
        printCurrents();
        printVoltmeters();
    }
};

void program() {
    int choice;
    cout << "\n--------------------------------------------------------------\n"
            "DCCalculator - Aplikacija za proračun DC mreže korištenjem MNČ\n"
            "--------------------------------------------------------------\n\n"
            "1 - Unos preko programa (grana po grani sa vrijednostima)\n"
            "2 - Unos preko Falstad simulatora  (nakon crtanje mreže i\n"
            "    postavljanje komponenti, Export As Text...)\n\n"
            "Vaš izbor: ";
    cin >> choice;
    cout << endl;
    if (choice == 2) {
        ifstream f("falstad.txt");
        if (!f.good()) {
            ofstream file("falstad.txt");
            file.close();
        }
        cout << "Postavite mrežu u \"falstad.txt\" datoteku (datoteka je kreirana\n"
                "ako nije vec postojala), spasite izmjene u datoteku i kliknite\n"
                "enter za nastavak...";
        cin.ignore();
        cin.ignore();
        cout << endl;
        Circuit cir("falstad.txt");
        cir.setRefNode(1);
        cir.solve();
        cir.printSolution();
        return;
    }
    ofstream file("falstad.txt");
    cout << "Prvo se unese broj čvorova (dijelovi  mreže  koji  nemaju isti\n"
            "potencijal,  tretirati   ampermetar  i  vatmetar   kao  obicne\n "
            "komponente), onda  se unese broj grana između svaka dva čvora,\n"
            "nakon toga  se unese šta svaka grana sadrži i vrijednost toga.\n"
            "Kliknite enter za nastavak...";
    cin.ignore();
    cin.ignore();
    int n, m, iTemp, jTemp, kTemp, type, refNode, vNode;
    double value;
    vector<Branch> branches;
    cout << "\nBroj čvorova u mreži: ";
    cin >> n;
    for (int i(1); i <= n - 1; i++) {
        for (int j(i + 1); j <= n; j++) {
            cout << "\nBroj grana između " << i << ". i " << j << ". čvor (0 ako nisu povezani): ";
            cin >> m;
            for (int k(1); k <= m; k++) {
                if(m == 1)
                    cout << "Grana između " << i << ". i " << j << ". čvor sadrži (R - 1, E - 2, Is - 3, Uv - 4, Ia - 5, Pw - 6): ";
                else
                    cout<< k << ". grana između " << i << ". i " << j << ". čvor sadrži (R - 1, E - 2, Is - 3, Uv - 4, Ia - 5, Pw - 6): ";
                cin >> type;
                if      (type == 1) cout << "Vrijednosti otpornika: ";
                else if (type == 2) cout << "Vrijednosti naponskog izvora (pomnožiti sa -1 ako je + na strani čvora " << i <<"): ";
                else if (type == 3) cout << "Vrijednosti strujnog izvora (pomnožiti sa -1 ako strelica ide ka čvoru " << i <<"): ";
                else if (type == 4) cout << "Ako je plus strana voltmetra na strani čvora " << j << " unesite -1, inače 1: ";
                else if (type == 5) cout << "Ako strelica ide ka čvoru " << i << " unesite -1, inače 1: ";
                else if (type == 6) {
                    cout << "Ako je plus strana vatmetra (za struju) na strani čvora " << j << " (obrnut) unesite -1, inače 1: ";
                    cin >> value;
                    iTemp = i;
                    jTemp = j;
                    if (value == -1) {
                        jTemp = i;
                        iTemp = j;
                    }
                    Branch branch(iTemp, jTemp, 0, 8, 0);
                    appendToFile("falstad.txt", branch);
                    cout << "Sa kojim trećim čvorom je vatmetar povezan (za mjerenje napona): ";
                    again:
                    cin >> vNode;
                    if (vNode > n) {
                        cout << "Broj čvora je izvan granice, unesite opet: ";
                        goto again;
                    }
                    cout << "Ako je plus strana vatmetra (za napona) na strani čvora " << vNode << " (obrnut) unesite -1, inače 1: ";
                }
                cin >> value;
                iTemp = i;
                jTemp = j;
                if (type == 6) {
                    type++;
                    jTemp = vNode;
                }
                if ((type == 4 || type == 5 || type == 6) && value == -1) {
                    int temp = iTemp;
                    iTemp = jTemp;
                    jTemp = temp;
                }
                if (m == 1) kTemp = 0;
                else kTemp = k;
                Branch branch(iTemp, jTemp, kTemp, type, value);
                appendToFile("falstad.txt", branch);
            }
        }
    }
//    cout << "\nReferentni čvor: ";
//    cin >> refNode;
    refNode = 1;
    Circuit cir("falstad.txt");
    cir.setRefNode(refNode);
    cir.solve();
    cir.printSolution();
}

int main() {
    try {
        program();
//  For faster testing comment line above and uncomment lines below
//        Circuit cir("falstad.txt");
//        cir.setRefNode(1);
//        cir.solve();
//        cir.printSolution();
    } catch (exception &e) {
        cerr << e.what();
        program();
    }
}
