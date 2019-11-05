#include "Topology.hpp"

using std::cout;
using std::endl;

int
main()
{
  Topology t1(6), t2(6);

  t1.insert_edge(0,1);
  t1.insert_edge(2,3);
  t1.insert_edge(4,5);
  t1.insert_edge(1,3);
  t1.insert_edge(3,5);
  t1.insert_edge(5,0);
  t1.insert_edge(0,2);
  t1.insert_edge(2,4);
  t1.insert_edge(4,1);

  t2.insert_edge(0,4);
  t2.insert_edge(4,5);
  t2.insert_edge(5,2);
  t2.insert_edge(2,1);
  t2.insert_edge(1,3);
  t2.insert_edge(3,0);
  t2.insert_edge(1,4);
  t2.insert_edge(0,2);
  t2.insert_edge(3,5);

  cout << t1.node_labelling();
  cout << t2.node_labelling();

  if (isomorphic(t1,t2))
    cout << "isomorphic" << endl;
  else
    cout << "not isomorphic" << endl;

  vector<int> perm = t2.node_labelling();

  for (int i = 0; i < t2.n_nodes(); ++i)
    {
      for (int j = 0; j < t2.n_nodes(); ++j)
	cout << t2.adjacency(perm[i],perm[j]) << " ";
      cout << endl;
    }
}
