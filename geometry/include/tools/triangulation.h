#pragma once
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unordered_set>
#include "declare.h"
#include "kdtree.h"

namespace tnurbs
{
	struct  VertexInsertionOrder
	{
		/**
		 * The Enum itself
		 * @note needed to pre c++11 compilers that don't support 'class enum'
		 */
		enum Enum
		{
			/**
			 * Automatic insertion order optimized for better performance
			 * @details breadth-first traversal of a Kd-tree for initial bulk-load,
			 * randomized for subsequent insertions
			 */
			Auto,
			/// insert vertices in same order they are provided
			AsProvided,
		};
	};


	template <typename T, int dim>
	struct OpenMesh::vector_traits<Eigen::Vector<T, dim>>
	{
		using vector_type = Eigen::Vector<T, dim>;
		using value_type = T;
		static const size_t size_ = dim;
		static size_t size() { return size_; }
	};




	template<typename T, int dim>
	struct TriTraits : public OpenMesh::DefaultTraits
	{
		typedef Eigen::Vector<T, dim> Point;
		// typedef Eigen::Vector<T, dim> Point;
	};

	template<typename T, int dim>
	using Mesh = OpenMesh::TriMesh_ArrayKernelT<>;
	// using Mesh = OpenMesh::PolyMesh_ArrayKernelT<TriTraits<T, dim>>;
	
	//现在dim必须为2，如果将kdtree拓展到任意维空间，既可以将dim的限制放开
	template<typename T, int dim = 2, typename TNearPointLocator = KDTree<T>>
	class triangulation
	{
		// typedef OpenMesh::PolyMesh_ArrayKernelT<TriTraits<T, dim>>::EdgeHandle EdgeHandle;
		typedef Mesh<T, dim>::EdgeHandle EdgeHandle ;
		typedef TriTraits<T, dim>::Point Point;
		typedef Point Vector;
		typedef Mesh<T, dim>::FaceHandle FaceHandle;
		typedef Mesh<T, dim>::VertexHandle VertexHandle;
		typedef Mesh<T, dim>::HalfedgeHandle HalfedgeHandle;
		typedef Mesh<T, dim>::FaceHalfedgeIter FaceHalfedgeIter;
		typedef Mesh<T, dim>::VertexOHalfedgeIter VertexOHalfedgeIter;
			
		struct PositionType
		{
			/**
			 * The Enum itself
			 * @note needed to pre c++11 compilers that don't support 'class enum'
			 */
			enum Enum
			{
				InFace,
				OnEdge,
				OnVertex,
			};
		};


		struct PointOnFace
		{
			PositionType pos_type;
			VertexHandle vh;
			EdgeHandle eh;
			FaceHandle fh;
		};
	private:

		Mesh<T, dim> m_mesh;
		std::vector<Point> m_fixed_edges;
		std::vector<Point> m_vertices;

		VertexInsertionOrder m_vertex_insertion_order;
		TNearPointLocator m_nearPt_locator;

		// Eigen::Vector<T, dim>
	public:
		triangulation() = default;
		bool is_finalized()
		{
			// unimplement
			return false;
		}

		void add_super_triangle(const Box<T, dim>& box)
		{
			// unimplement
			return;
		}


		void add_super_triangle(Box<T, dim>& box)
		{
			// unimplement
			return;
		}

		PointOnFace walk_triangles(const VertexHandle& start_vertex, const Point& position) const
		{
			bool found = false;
			PointOnFace result;
			HalfedgeHandle half_edge = m_mesh.halfedge_handle(start_vertex);
			while (found == false)
			{
				found = true;
				for (FaceHalfedgeIter it = m_mesh.cfh_iter(half_edge); it.is_valid(); ++it)
				{
					HalfedgeHandle pre_halfedge = m_mesh.prev_halfedge_handle(it.current_halfedge_handle);
					VertexHandle pre_veterx = m_mesh.from_vertex_handle(pre_halfedge);
					VertexHandle current_vertex = m_mesh.from_vertex_handle(it.current_halfedge_handle());
					VertexHandle next_vertex = m_mesh.to_vertex_handle(it.current_halfedge_handle());

					Point pre_point = m_mesh.point(pre_veterx);
					Point current_point = m_mesh.point(current_vertex);
					Point next_point = m_mesh.point(next_vertex);
				
					Vector e1 = pre_point - current_point;
					Vector e2 = next_point - current_point;

					Vector e = position - current_point;

					if (e.norm() < PRECISION<T>::value)
					{
						// on vertex
						result.pos_type = PositionType::OnVertex;
						result.vh = current_vertex;
						return result;
					}
					if ((position - next_point).norm() < PRECISION<T>::value)
					{
						// on vertex
						result.pos_type = PositionType::OnVertex;
						result.vh = next_vertex;
						return result;
					}
					if (e.dot(e2) / e2.norm() < PRECISION<T>::value)
					{
						if (e.norm() < e2.norm() && e.dot(e2) > 0)
						{
							//on edge
							result.pos_type = PositionType::OnEdge;
							result.eh = m_mesh.edge_handle(half_edge);
							return result;
						}

						// point is on curve of edge, but not on edge
						continue;
					}

					Eigen::Vector2<T> vec;
					vec[0] = e.dot(e1);
					vec[1] = e.dot(e2);
					Eigen::Matrix<T, dim, 2> mat;
					mat(0, 0) = e1.dot(e1);
					mat(0, 1) = mat(1, 0) = e1.dot(e1);
					mat(1, 1) = e1.dot(e1);
					Eigen::JacobiSVD<Eigen::Matrix<T, dim, 2>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
					Eigen::Vector<T, 2> coeffs = matSvd.solve(vec);

					if (matSvd.info() != Eigen::Success)
					{
						assert(false);
						return ENUM_NURBS::NURBS_ERROR;
					}
					if (coeffs[0] < 0)
					{
						// 点在edge的右侧(即：外侧)
						if (m_mesh.is_boundary(half_edge) == true)
						{
							assert(false);
						}
						half_edge = m_mesh.opposite_halfedge_handle(it.current_halfedge_handle());
						found = false;
						break;
					}
				}
			}
			result.pos_type = PositionType::InFace;
			result.fh = m_mesh.face_handle(half_edge);
			return result;
		}

		VertexHandle insert_vertex_inside_triangle(Point point, FaceHandle fh)
		{
			VertexHandle vh = m_mesh.add_vertex(point);
			m_mesh.split(fh, vh);
			return vh;
		}
		VertexHandle insert_vertex_on_edge(Point point, EdgeHandle eh)
		{
			VertexHandle vh = m_mesh.add_vertex(point);
			m_mesh.split_edge(eh, vh);
			return vh;
		}
		void ensure_Delaunay_by_edge_flips(VertexHandle vh)
		{
			VertexOHalfedgeIter it = m_mesh.voh_iter(vh);
			for (; it.is_valid(); ++it)
			{
				HalfedgeHandle hh = it.current_halfedge_handle();
				bool is_bound = m_mesh.is_boundary(hh);
				assert(is_bound == false);
				HalfedgeHandle next_hh = m_mesh.next_halfedge_handle(hh);
				std::vector<HalfedgeHandle> flip_edges{ next_hh };
				while (flip_edges.empty())
				{
					HalfedgeHandle current_hh = flip_edges.back();
					flip_edges.pop_back();
					bool can_flip = m_mesh.is_flip_ok(m_mesh.edge_handle(current_hh));
					assert(can_flip == true);
					m_mesh.flip(m_mesh.edge_handle(current_hh));
				}
			}
			return;
		}
		
		void insert_vertex(const Point& position, const VertexHandle& walk_start_vertex)
		{
			PointOnFace ptf = walk_triangles(walk_start_vertex, position);
			VertexHandle vh;
			switch (ptf.pos_type)
			{
			case PositionType::OnVertex:
				assert(false);
				break;
			case PositionType::OnEdge:
				vh = insert_vertex_on_edge(position, ptf.eh);
				break;
			case PositionType::InFace:
				vh = insert_vertex_inside_triangle(position, ptf.fh);
				break;
			default:
				break;
			}
			ensure_Delaunay_by_edge_flips(vh);
			return;
		}
		
		
		void insert_vertices_KDTreeBFS(size_t exits_vertex_count, size_t new_vertex_count,const Box<T, dim>& box, const std::vector<Point>& vertices)
		{
			// unimplement
			for (size_t index = 0; index < new_vertex_count; ++index)
			{
				// VertexHandle vertex_handle = m_mesh.add_vertex(vertex);
				// vertex_handle.idx();

			}
			return;
		}
		void insert_vertices_randomized()
		{
			// unimplement
			return;
		}
		bool insert_vertices(const std::vector<Point>& vertices)
		{
			// if (is_finalized())
			// {
			// 	return false;
			// }
			bool is_first_time = m_mesh.vertices_empty();
			const T max = std::numeric_limits<T>::max();
			Box<T, dim> box;
			box.Min.setConstant(-max);
			box.Max.setConstant(max);
			if (is_first_time)
			{
				box = envelop_points(vertices);
				add_super_triangle(box);
			}
			{
				if (m_nearPt_locator.empty() == true && m_vertices.empty() == false)
				{
					Box<T, dim> box;
					box.Min.setConstant(-max);
					box.Max.setConstant(max);
					box = envelop_points(m_vertices);
					m_nearPt_locator = KDTree<T>(box);
					for (const auto& point : m_vertices)
					{
						m_nearPt_locator.insert(point);
					}
				}
			}
			size_t exits_vertex_count = m_vertices.size();
			size_t new_vertex_count = vertices.size();
			size_t vertices_count = exits_vertex_count + new_vertex_count;
			for (const auto& vertex : vertices)
			{
				m_vertices.push_back(vertex);
			}
			switch (m_vertex_insertion_order)
			{
			case VertexInsertionOrder::AsProvided:
				//TODO:
				assert(false);
				break;
			case VertexInsertionOrder::Auto:
				is_first_time ? insert_vertices_KDTreeBFS(exits_vertex_count, new_vertex_count, box, vertices)
					: insert_vertices_randomized();
			default:
				break;
			}
			return true;
		}




	};
}


