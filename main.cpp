#include <iostream>
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using namespace std;

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
    vector<Branch> branches;

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
                e(k, 0) = b.getValue();
                k++;
            }
        }
        return e;
    }

public:
    Circuit(int noOfNodes, int refNode) : noOfNodes(noOfNodes), refNode(refNode) { solved = 0; }

    int getNoOfNodes() const { return noOfNodes; }
    int getRefNode() const { return refNode; }
    const vector<Branch> &getBranches() const { return branches; }
    const vector<double> &getNodesVoltages() const { return nodesVoltages; }
    const vector<double> &getVolSourcesCurrents() const { return volSourcesCurrents; }

    void setNoOfNodes(int noOfNodes) { Circuit::noOfNodes = noOfNodes; solved = 0; }
    void setRefNode(int refNode) { Circuit::refNode = refNode; solved = 0; }
    void setBranches(const vector<Branch> &branches) { Circuit::branches = branches; solved = 0; }

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

        solved = 1;
    }
};

int main() {
}
