/*
 * disjoint_set.hpp
 * Author: Aven Bross
 *
 * Implementation of a simple disjoint set data structure.
 */

#ifndef __DISJOINT_SET_HPP
#define __DISJOINT_SET_HPP

// STL headers
#include <vector>

/*
 * disjoint_set
 * This structure provides a simple disjoint set structure with union and find operations.
 * New singleton sets are constructed in order using the make_next() function. Union,
 * find, comparison, and existance functions are provided.
 */

class disjoint_set {
	public:
		disjoint_set() : next(0) {}
		
		int make_next() {
			sets.push_back(set_data(next));
			return next++;
		}
		
		bool exists(int n) {
			return  (n >= 0) && (std::size_t) n < sets.size();
		}
		
		int find(int n) {
			if(sets[n].parent != n) {
				sets[n].parent = find(sets[n].parent);
			}
			return sets[n].parent;
		}
		
		bool compare(int n, int m) {
			if(!exists(n) || !exists(m)) return false;
			return find(n) == find(m);
		}
				
		int take_union(int n, int m) {
			if(!exists(n)) return m;
			if(!exists(m)) return n;
			
			int n_root = find(n), m_root = find(m);
			
			if(n_root == m_root) return n_root;
			
			if(sets[n_root].rank < sets[m_root].rank) {
				sets[n_root].parent = m_root;
				return m_root;
			}
			else if(sets[n_root].rank > sets[m_root].rank) {
				sets[m_root].parent = n_root;
				return n_root;
			}
			else {
				sets[m_root].parent = n_root;
				++sets[n_root].rank;
				return n_root;
			}
		}
		
	private:
		struct set_data {
			int parent, rank;
			set_data(int parent) : parent(parent), rank(0) {}
		};
		
		std::vector<set_data> sets;
		int next;
};

#endif