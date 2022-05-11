#include <iostream>
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using namespace std;

class Branch {
    int nodeI, nodeJ, branchK, type; // type = 1 for R, type = 2 for E
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
    int noOfNodes, refNode;
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
    MatrixXd admittances() { return MatrixXd::Zero(1, 1); }
    MatrixXd volSourcesConnections() { return MatrixXd::Zero(1, 1); }
    MatrixXd currentSources() { return MatrixXd::Zero(1, 1); }
    MatrixXd voltageSources() { return MatrixXd::Zero(1, 1); }

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
        int n = noOfNodes - 1, m = noOfVolSources();
        MatrixXd G(n, n), B(n, m), D = MatrixXd::Zero(m, m), A(n+m, n+m), i(n, 1), e(m, 1), z(n+m, 1), x;
        G = admittances();
        B = volSourcesConnections();
        A << G,             B,
             B.transpose(), D;
        i = currentSources();
        e = voltageSources();
        z << i, e;
        x = A.inverse() * z;

        solved = 1;
    }
};

int main() {
}
