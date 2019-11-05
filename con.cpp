#include <iterator>
#include "Topology.hpp"

using std::ostream_iterator;
using std::cout;
using std::endl;

int
main()
{
  Topology t(9);

  t.insert_edge(5,1);
  t.insert_edge(1,4);
  t.insert_edge(1,0);
  t.insert_edge(4,0);
  t.insert_edge(4,3);
  t.insert_edge(0,2);
  t.insert_edge(2,3);
  t.insert_edge(3,6);
  t.insert_edge(2,7);
  t.insert_edge(2,8);
  t.insert_edge(7,8);
  t.insert_edge(7,8);

  t.postscript_print("con.ps");

  vector<TopologyComponent> components = t.biconnected_components();

  for (vector<TopologyComponent>::iterator c = components.begin();
       c != components.end(); ++c)
    {
      cout << "nodes: ";
      copy(c->_nodes.begin(), c->_nodes.end(),
	   ostream_iterator<int>(cout, " "));
      if (c->_vacuum) cout << ", vacuum";
      cout << endl;
    }
}
