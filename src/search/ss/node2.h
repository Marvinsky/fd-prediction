#ifndef NODE2_H
#define NODE2_H

#include <vector>

class Node2 {
private:
      int f_value;
      int level;
      int q_value;
public:
      Node2();
      Node2(int f, int l);
      int getF();
      void setF(int f);
      int getL();
      void setL(int l);
      int getQ();
      void setQ(int q);
      friend bool operator< (const Node2 &n1, const Node2 &n2);
};
#endif
