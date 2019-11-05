/******************************************************************************
 *                                                                            *
 * Copyright (C) 2002-2005 Michal Czakon                                      *
 *                                                                            *
 ******************************************************************************/

#include "TopologyGenerator.hpp"

using namespace std;

/******************************************************************************
 *                                                                            *
 * next_topology                                                              *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * If the first user call after definition of a TopologyGenerator object and
 * setting of options is to next_topology(), then the method returns true if
 * there are any topologies consistent with settings or false otherwise. In the
 * current implementation, this method either moves the pointer within the list
 * of topologies or makes a new one when necessary.
 *
 */

bool
TopologyGenerator::next_topology()
{
  // initialize if it's the first call

  if (_topologies.empty())
    {
      while (!make_topologies())
	if (_fixed_node_count || _node_count[2] < 2)
	  return false;
	else
	  {
	    _node_count[2] -= 2;
	    ++_node_count[3];
	  }	      

      _current_topology = _topologies.begin();
      return true;
    }

  if (++_current_topology != _topologies.end()) return true;

  // go to the next node partition if at the end of topology list
  
  if (_fixed_node_count)
    {
      _current_topology = _topologies.begin();
      return false;
    }
  else
    {
      _topologies.clear();

      do
	{
	  if (_node_count[2] < 2)
	    {
	      _node_count[2] += 2*_node_count[3];
	      return false;
	    }
	  else
	    {
	      _node_count[2] -= 2;
	      ++_node_count[3];
	    }
	}
      while (!make_topologies());
	  
      _current_topology = _topologies.begin();
      return true;
    }
}

/******************************************************************************
 *                                                                            *
 * make_topologies                                                            *
 *                                                                            *
 ******************************************************************************/

bool
TopologyGenerator::make_topologies()
{
  first_adjacency();

  do
    {
      Topology t(_n_nodes);

      for (int parent_node = _n_nodes-1; parent_node >= 0; --parent_node)
	{
	  for (int child_node = 0; child_node < parent_node; ++child_node)
	    for (int count = _adjacency[parent_node][child_node];
		 count > 0; --count)
	      {
		t.insert_edge(child_node, parent_node);
	      }

	  for (int count = _unassigned[parent_node][parent_node]/2;
	       count > 0; --count)
	    {
	      t.insert_edge(parent_node, parent_node);
	    }
	}

      if ((_options & AllowDisconnected || t.is_connected())
	  && 
	  !(_options& OneParticleIrreducible&& !t.is_one_particle_irreducible())
	  &&
	  !(_options & NoTadpoles && t.has_tadpoles())
	  &&
	  !(_options & NoSelfEnergies && t.has_self_energies())
	  &&
	  !(_options & OnShell && !t.is_on_shell()))
	{
	  if (!(_options & EquivalentExternalNodes)) t.fix_external_nodes();
	  _topologies.insert(t);
	}
    }
  while (next_adjacency());

  return !_topologies.empty();
}

/******************************************************************************
 *                                                                            *
 * first_adjacency                                                            *
 *                                                                            *
 ******************************************************************************/

void
TopologyGenerator::first_adjacency()
{
  vector<int> unassigned;
  for (int i = 0; i < _node_count.size(); ++i)
    unassigned.insert(unassigned.end(), _node_count[i], i+1);
  
  _n_nodes = unassigned.size();
  _sort = vector<vector<bool> >(_n_nodes, vector<bool>(_n_nodes, true));
  _compare = vector<vector<bool > >(_n_nodes, vector<bool>(_n_nodes, false));
  _unassigned = vector<vector<int> >(_n_nodes, vector<int>(_n_nodes));
  _adjacency = vector<vector<int> >(_n_nodes, vector<int>(_n_nodes));

  _unassigned[_n_nodes-1] = unassigned;
  
  int node = 0;
  vector<int>::const_iterator count = _node_count.begin();

  if (_options & EquivalentExternalNodes)
    {
      _sort[_n_nodes-1][0] = false;
       if (node += *count) _compare[node-1].clear();
    }
  else
    for (; node < *count; ++node)
      {
	_sort[_n_nodes-1][node] = false;
	_compare[node].clear();
      }

  for (++count; node < _n_nodes; ++count)
    {
      _sort[_n_nodes-1][node] = false;
      if (node += *count) _compare[node-1].clear();
    }

  for (int node = _n_nodes-1; node > 0; --node)
    first_connection(node);
}

/******************************************************************************
 *                                                                            *
 * next_adjacency                                                             *
 *                                                                            *
 ******************************************************************************/

bool
TopologyGenerator::next_adjacency()
{
  int node;
  do
    {
      for (node = 1; node < _n_nodes && !next_connection(node); ++node) ;
      if (node == _n_nodes) return false;
      for (--node; node > 0 && first_connection(node); --node) ;
    }
  while (node || (_unassigned[0][0] & 1) || !is_valid(0));

  return true;
}

/******************************************************************************
 *                                                                            *
 * first_connection                                                           *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Connects a node to its children (nodes with indices smaller than the current)
 * If the connection is not valid, then calls next_connection(). If successful
 * recalculates the validation graph, i.e. sets _compare for all nodes which
 * should be compared with the current.
 *
 */

bool
TopologyGenerator::first_connection(int node)
{
  connect(node, 0);
  
  if ((_unassigned[node][node] & 1) || !is_valid(node))
    return next_connection(node);

  // recalculate the validation graph

  bool compare = _sort[node][node];
  for (int child_node = node-1; child_node >= 0 &&
	 !_compare[child_node].empty(); --child_node)
    {
      _compare[child_node][node] = compare;
      if (!_sort[node][child_node]) compare = false;
    }

  return true;
}

/******************************************************************************
 *                                                                            *
 * next_connection                                                            *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * The algorithm generates in fact succeding partitions of the remaining node
 * degree among the children. However the partitions are sorted with respect
 * to identical children nodes and the node is validated according to the
 * validation graph (see is_valid()). If a valid node is found the validation
 * graph is recalculated.
 *
 */

bool
TopologyGenerator::next_connection(int node)
{
  do
    {
      int child_node;
      for (child_node = node-1; child_node >= 0 &&
	     !_adjacency[node][child_node]; --child_node) ;
      
      if (child_node == -1) return false;
      
      --_adjacency[node][child_node];

      ++_unassigned[node][node];
      ++_unassigned[node][child_node];

      ++_unassigned[node-1][child_node];
      _sort[node-1][child_node] = false;
      
      connect(node, child_node+1);
    }
  while ((_unassigned[node][node] & 1) || !is_valid(node));

  // recalculate the validation graph

  bool compare = _sort[node][node];
  for (int child_node = node-1; child_node >= 0 &&
	 !_compare[child_node].empty(); --child_node)
    {
      _compare[child_node][node] = compare;
      if (!_sort[node][child_node]) compare = false;
    }

  return true;
}

/******************************************************************************
 *                                                                            *
 * connect                                                                    *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Connects the remaining legs of a node to children with indices starting from
 * first_node. The connections are sorted among identical nodes. Updates also
 * _sort and _unassigned of the node's child.
 *
 */

void
TopologyGenerator::connect(int node, int first_node)
{
  for (int child_node = first_node; child_node < node; ++child_node)
    {
      int value =
	min(_unassigned[node][child_node], _unassigned[node][node]);

      if (_sort[node][child_node] && value >= _adjacency[node][child_node-1])
	{
	  value = _adjacency[node][child_node-1];
	  _sort[node-1][child_node] = true;
	}
      else
	{
	  _sort[node-1][child_node] = false;
	}
      
      _adjacency[node][child_node] = value;

      _unassigned[node][node] -= value;
      _unassigned[node][child_node] -= value;

      _unassigned[node-1][child_node] = _unassigned[node][child_node];
    }
}

/******************************************************************************
 *                                                                            *
 * is_valid                                                                   *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Checks whether a node is valid according to the validation graph, i.e.
 * whether it is not greater than its parent in the graph. As a speedup in 
 * generation of connected topologies a test is also performed whether external
 * nodes have to be connected.
 *
 */

bool
TopologyGenerator::is_valid(int node) const
{
  if (!(_options & AllowDisconnected) && node == _node_count[0])
    for (int child_node = 0; child_node < node; ++child_node)
      if (_unassigned[node][child_node]) return false;
  
  if (_compare[node].empty()) return true;
  
  for (int parent_node = node+1; parent_node < _n_nodes; ++parent_node)
    if (_compare[node][parent_node] && compare(parent_node, node) < 0)
      return false;

  return true;
}

/******************************************************************************
 *                                                                            *
 * compare                                                                    *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * parent_node is compared with child_node lexicographically after sorting the
 * child node adjacencies within symmetry cells of the parent. The
 * implementation of sorting is a version of the bucket sort for speed
 * (although it was not proved that using multiset would hamper the overall
 * performance)
 *
 */

int
TopologyGenerator::compare(int parent_node, int child_node) const
{
  int node = 0;
  do
    {
      const int max_adjacency = _node_count.size()+1;
      vector<int> parent_bucket(max_adjacency);
      vector<int> child_bucket(max_adjacency);
      int max_parent_adjacency = 0;
      int max_child_adjacency = 0;

      do
	{
	  const int parent_adjacency = _adjacency[parent_node][node];

	  if (node != parent_node)
	    {
	      ++parent_bucket[parent_adjacency];
	      if (parent_adjacency > max_parent_adjacency)
		max_parent_adjacency = parent_adjacency;
	    }

	  int child_adjacency = _adjacency[child_node][node];
	  if (node > child_node) child_adjacency = _adjacency[node][child_node];
	  
	  if (node != child_node)
	    {
	      ++child_bucket[child_adjacency];
	      if (child_adjacency > max_child_adjacency)
		max_child_adjacency = child_adjacency;
	    }
	}
      while (node < parent_node && _sort[parent_node][++node]);

      if (max_parent_adjacency != max_child_adjacency)
	return max_parent_adjacency-max_child_adjacency;
      
      for (int i = max_parent_adjacency; i > 0; --i)
	if (parent_bucket[i] != child_bucket[i])
	  return parent_bucket[i] - child_bucket[i];
    }
  while (node < parent_node);

  return 0;
}
