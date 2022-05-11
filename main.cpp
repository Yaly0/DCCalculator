#include <iostream>
#include <Eigen/Dense>

using namespace std;

class Branch {
    int nodeI, nodeJ, branchK, type; // type = 1 for R, type = 2 for E
    double value;
public:
    Branch(int nodeI, int nodeJ, int branchK, int type, double value) : nodeI(nodeI), nodeJ(nodeJ), branchK(branchK),
                                                                        type(type), value(value) {}
    int getNodeI() const { return nodeI; }
    void setNodeI(int nodeI) { Branch::nodeI = nodeI; }
    int getNodeJ() const { return nodeJ; }
    void setNodeJ(int nodeJ) { Branch::nodeJ = nodeJ; }
    int getBranchK() const { return branchK; }
    void setBranchK(int branchK) { Branch::branchK = branchK; }
    int getType() const { return type; }
    void setType(int type) { Branch::type = type; }
    double getValue() const { return value; }
    void setValue(double value) { Branch::value = value; }
};

int main() {
}
