#include <iostream>
#include <vector>
#include <Eigen/Dense>

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

        solved = 1;
    }
};

int main() {
}
