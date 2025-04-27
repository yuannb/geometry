
// �ο�(cao) https://mozilla.org/MPL/2.0/.
#pragma once
#include <limits>
#include <array>
#include <vector>
#include "declare.h"
#include <assert.h>
namespace tnurbs
{
	struct node_split_direction
	{
		enum Enum
		{
			X,
			Y,
		};
	};

	template<typename number_type, size_t number_vertices_in_leaf = 32>
	class KDTree
	{
	public:
		typedef unsigned int node_index;
		typedef unsigned int point_index;
		typedef Eigen::Vector2<number_type> point_type;
		typedef std::array<int, 2> children_type;
	private:
		struct Node
		{
			children_type children; ///< two children if not leaf; {0,0} if leaf
			std::vector<point_index> data; ///< points' data if leaf
			Node() 
			{
				children[0] = 0;
				children[1] = 0;
			}
			void set_children(const node_index c1, const node_index c2)
			{
				children[0] = 0;
				children[1] = 0;
			}
			bool is_leaf() const
			{
				return children[0] == children[1];
			}
		};

		struct nearest_task
		{
			node_index node;
			node_split_direction::Enum dir;
			number_type dist2;
			Box<number_type, 2> box;
			nearest_task() = default;
			nearest_task(const node_index node_,
				const node_split_direction::Enum dir_,
				const number_type dist2_,
				const Box<number_type, 2> box_)
				: node(node_)
				, dir(dir_)
				, dist2(dist2_)
				, box(box_)
			{}
		};

	private:
		node_index m_root;
		std::vector<Node> m_nodes;
		Box<number_type, 2> m_box;
		node_split_direction::Enum m_root_dir;
		number_type m_tolerance;
		bool m_is_root_box_initialized;
		
		mutable std::vector<nearest_task> m_tasks_stack;
		
		node_index add_new_node()
		{
			const node_index new_node_index = static_cast<node_index>(m_nodes.size());
			m_nodes.push_back(Node());
			return new_node_index;
		}

		std::size_t which_child(const point_type& point, const number_type split, const node_split_direction::Enum dir) const
		{
			return static_cast<std::size_t>(
				dir == node_split_direction::X ? point[0] > split : point[1] > split);
		}

		static void calc_split_info(const Box<number_type, 2>& box, node_split_direction::Enum dir, number_type& mid
			, node_split_direction::Enum& new_dir)
		{
			point_type point = box.get_middle_point();
			switch (dir)
			{
			case tnurbs::node_split_direction::X:
				mid = point[0];
				return;
			case tnurbs::node_split_direction::Y:
				mid = point[1];
				return;
			}
			return;
		}

		//extendһ�β�һ��Ϊ����point
		bool extend_tree(const point_type& point)
		{
			if (m_box.is_contain_point(point, m_tolerance) == true)
			{
				assert(false);
				return false;
			}
			const node_index new_root = add_new_node();
			const node_index new_leaf = add_new_node();
			switch (m_root_dir)
			{
			case tnurbs::node_split_direction::X:
				m_root_dir = node_split_direction::Y;
				if (point[1] < m_box.Min[1])
				{
					m_nodes[new_root].set_children(new_leaf, m_root);
					m_box.Min[1] -= (m_box.Max[1] - m_box.Min[1]);
				}
				else if (point[1] > m_box.Max[1])
				{
					m_nodes[new_root].set_children(m_root, new_leaf);
					m_box.Max[1] += (m_box.Max[1] - m_box.Min[1]);
				}
				else
				{
					assert(false);
				}
				break;
			case tnurbs::node_split_direction::Y:
				m_root_dir = node_split_direction::X;
				if (point[0] < m_box.Min[0])
				{
					m_nodes[new_root].set_children(new_leaf, m_root);
					m_box.Min[0] -= (m_box.Max[0] - m_box.Min[0]);
				}
				else if (point[0] > m_box.Max[0])
				{
					m_nodes[new_root].set_children(m_root, new_leaf);
					m_box.Max[0] += (m_box.Max[0] - m_box.Min[0]);
				}
				else
				{
					assert(false);
				}
				break;
			}
			m_root = new_root;
			return true;
		}

		//ֻ����kdtreeֻ��һ���ڵ��ʱ�����
		void initialize_root_box(const std::vector<point_type>& points)
		{
			const std::vector<point_index>& data = m_nodes[m_root].data;
			m_box = Box<number_type, 2>(points[data.front()], points[data.front()]);
			for (auto it = data.begin() + 1; it != data.end(); ++it)
			{
				m_box.enlarge(points[*it]);
			}

			// Make sure bounding box does not have a zero size by adding padding:
			// zero-size bounding box cannot be extended properly
			if (m_box.Max[0] - m_box.Min[0] < m_tolerance)
			{
				m_box.Max[0] += 1;
			}
			if (m_box.Max[1] - m_box.Min[1] < m_tolerance)
			{
				m_box.Max[1] += 1;
			}
			m_is_root_box_initialized = true;
			return;
		}
	public:
		KDTree()
			: m_root_dir(node_split_direction::X)
			, m_box(Eigen::Vector2<number_type>(std::numeric_limits<number_type>::max(), std::numeric_limits<number_type>::max())
				, Eigen::Vector2<number_type>(std::numeric_limits<number_type>::max(), std::numeric_limits<number_type>::max()))
			, m_is_root_box_initialized(false)
		{
			m_root = add_new_node();
			m_tolerance = PRECISION<number_type>::value;
		};
		KDTree(const point_type& min, const point_type& max)
			: m_root_dir(node_split_direction::X)
			, m_box(min, max)
			, m_is_root_box_initialized(true)
		{
			m_root = add_new_node();
		}
		KDTree(const Box<number_type, 2>& box)
			: m_root_dir(node_split_direction::X)
			, m_box(box)
			, m_is_root_box_initialized(true)
		{
			m_root = add_new_node();
		}
		void insert(const point_index& ipoint, const std::vector<point_type>& points)
		{
			const point_type& pos = points[ipoint];
			while (m_box.is_contain_point(pos, m_tolerance))
			{
				extend_tree(pos);
			}
			node_index node = m_root;
			node_split_direction::Enum dir = m_root_dir;
			node_split_direction::Enum new_dir(node_split_direction::X);
			number_type mid{ 0 };
			Box<number_type, 2> new_box;
			Box<number_type, 2> box;
			while (true)
			{
				if (m_nodes[node].is_leaf() == true)
				{
					std::vector<point_index>& current_points = m_nodes[node].data;
					if (current_points.size() < number_vertices_in_leaf)
					{
						current_points.push_back(ipoint);
						return;
					}
					if (m_is_root_box_initialized == false)
					{
						initialize_root_box(points);
						box = m_box;
					}
					// split a full leaf node
					calc_split_info(box, dir, mid, new_dir);
					const node_index c1 = add_new_node(), c2 = add_new_node();
					Node& n = m_nodes[node];
					n.set_children(c1, c2);
					std::vector<point_index>& c1_data = m_nodes[c1].data;
					std::vector<point_index>& c2_data = m_nodes[c2].data;
					for (auto it = n.data.begin(); it != n.data.end(); ++it)
					{
						which_child(points[*it], mid, dir) == 0
							? c1_data.push_back(*it)
							: c2_data.push_back(*it);
					}
					n.data.clear();
				}
				else
				{
					calc_split_info(box, dir, mid, new_dir);
				}
				//add point to child
				const std::size_t ichild = which_child(points[ipoint], mid, dir);
				if (ichild == 0)
				{
					if (dir == node_split_direction::X)
					{
						box.Max[0] = mid;
					}
					else
					{
						box.Max[1] = mid;
					}
				}
				else
				{
					if (dir == node_split_direction::X)
					{
						box.Min[0] = mid;
					}
					else
					{
						box.Min[1] = mid;
					}
				}
				node = m_nodes[node].children[ichild];
				dir = new_dir;
			}
			return;
		}

		point_index nearest(const point_type& point, std::vector<point_type>& points) const
		{
			point_index out{ 0 };
			number_type min_dis2 = std::numeric_limits<number_type>::max();
			m_tasks_stack.emplace_back (m_root, m_root_dir, min_dis2, m_box);
			while (m_tasks_stack.empty() == false)
			{
				nearest_task t = std::move(m_tasks_stack.back());
				m_tasks_stack.pop_back();
				if (t.dist2 > min_dis2 + m_tolerance * m_tolerance)
				{
					continue;
				}
				const Node& n = m_nodes[t.node];
				if (n.is_leaf() == true)
				{
					for (auto it = n.data.begin(); it != n.data.end(); ++it)
					{
						const point_type& p = points[*it];
						number_type dist2 = (point - p).squaredNorm();
						if (dist2 < min_dis2)
						{
							min_dis2 = dist2;
							out = *it;
						}
					}
				}
				else
				{
					number_type mid{ 0 };
					node_split_direction::Enum new_dir;
					calc_split_info(t.box, t.dir, mid, new_dir);
					number_type dist_to_mid = t.dist2 == node_split_direction::X
						? (point[0] - mid)
						: (point[1] - mid);
					number_type dis_mid2 = dist_to_mid * dist_to_mid;
					std::size_t ichild = which_child(point, mid, t.dir);
					std::array<Box<number_type, 2>, 2> sub_box = t.box.split_at_middle(new_dir == 0 ? 0 : 1);
					if (ichild == 0)
					{
						m_tasks_stack.emplace_back(n.children[1], new_dir, dis_mid2, sub_box[1]);
						m_tasks_stack.emplace_back(n.children[0], new_dir, dis_mid2, sub_box[0]);
					}
					else
					{
						m_tasks_stack.emplace_back(n.children[0], new_dir, dis_mid2, sub_box[0]);
						m_tasks_stack.emplace_back(n.children[1], new_dir, dis_mid2, sub_box[1]);
					}
				}
			}
			return out;
		}
	
		bool empty() const
		{
			return !m_is_root_box_initialized;
		}
	};

	// void Initialize(const std::vector<point_type>& points)
	// {
	// 	point_type min = points.front();
	// 	number_type xmin = min[0];
	// 	number_type xmax = min[0];
	// 	number_type ymin = min[1];
	// 	number_type ymax = min[1];

	// 	for (const auto& point : points)
	// 	{
	// 		xmin = std::min(point[0], xmin);
	// 		xmax = std::max(point[0], xmax);
	// 		ymin = std::min(point[1], ymin);
	// 		ymax = std::max(point[1], ymax);
	// 	}
	// }
}

