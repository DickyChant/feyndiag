#include <sstream>
#include "TopologyGenerator.hpp"

using std::ostringstream;
using std::cout;
using std::endl;

int
main()
{
  int count = 0;

#if 1
  TopologyGenerator generator(2, 2, OneParticleIrreducible);
#else
  vector<int> node_count(4);
  node_count[0] = 0;
  node_count[1] = 0;
  node_count[2] = 4;
  node_count[3] = 0;

  TopologyGenerator generator(node_count,
			      EquivalentExternalNodes |
			      OneParticleIrreducible |
			      NoSelfEnergies);
#endif

  while (generator.next_topology())
    {
      ++count;

      Topology t = generator.current_topology();
      t.assign_momenta();
      t.print_edge_list();
      cout << endl;

      ostringstream name;
      name << "top" << count << ".ps";
      t.postscript_print(name.str());
    }
  cout << count << " generated topologies" << endl;
}
