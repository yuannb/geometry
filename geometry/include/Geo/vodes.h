#include "declare.h"
#include "nurbs_surface.h"
#include "interval_alogrithm.h"
#include <memory>
#include "bezier_curve_int.h"
#include <ctime>
#include <chrono>

// #define DSSI

namespace tnurbs
{
	inline static void printNurbsControlPoints(std::vector<Eigen::Matrix<double, 3, Eigen::Dynamic>>& points)
	{
		for (int index = 0; index < points.size(); ++index)
		{
			std::cout << "line*********************" << std::endl;
			std::cout << points[index] << std::endl;
		}
		std::cout << "end**************" << std::endl;
		return;
	}
	//对于大量计算点的位置，包围盒等。统一申请一些内存，减少内存重复申请和析构的耗时

	template<typename surface_type>
	class surf_compute
	{
		using T = typename surface_type::Type;
		static constexpr int dim = surface_type::dimension;
		static constexpr bool is_rational = surface_type::is_ratio;
		static constexpr int point_size = surface_type::point_size;
	public:
		surface_type* m_surface;
		Eigen::VectorX<T> auxiliary_u_knots_vector;
		int valid_auxiliary_u_knots_count;
		Eigen::VectorX<T> auxiliary_v_knots_vector;
		int valid_auxiliary_v_knots_count;
		std::vector<Eigen::Matrix<T, point_size, Eigen::Dynamic>> auxiliary_points;
		std::vector<Eigen::Matrix<T, point_size, Eigen::Dynamic>> auxiliary_points2;
		int m_u_degree;
		int m_v_degree;

		Eigen::VectorX<T> m_u_knots_vector;
		Eigen::VectorX<T> m_v_knots_vector;

		Eigen::MatrixX<T> u_ders_basis_funs_array;
		Eigen::MatrixX<T> v_ders_basis_funs_array;
		Eigen::VectorX<T> u_basis_funs_array;
		Eigen::VectorX<T> v_basis_funs_array;
		Eigen::MatrixX<T> ndu;
		Eigen::MatrixX<T> ndu_trans;
		Eigen::VectorX<T> left;
		Eigen::VectorX<T> right;
		Eigen::VectorX<T> iterArray;
		Eigen::ArrayX<T> coeff_denominator;;
		Eigen::Array<T, Eigen::Dynamic, 1> first_coeff;
		Eigen::Array<T, Eigen::Dynamic, 1> second_coeff;
		std::vector<T> gammas;
		std::array<Eigen::Array<T, Eigen::Dynamic, 1>, 2> arrays;
		Eigen::Matrix<T, point_size,  Eigen::Dynamic> temp;
		Eigen::Matrix<T, point_size, Eigen::Dynamic> RW;
		std::vector<Eigen::Matrix<T, point_size, Eigen::Dynamic>> RW2;
		surf_compute(surface_type* surf)
		{
			m_surface = surf;
			const auto& control_points = m_surface->get_control_points_ref();
			int v_count = control_points.size();
			int u_count = control_points[0].cols();
			m_u_degree = m_surface->get_u_degree();
			m_v_degree = m_surface->get_v_degree();
			auxiliary_points.resize(v_count + 2 * m_v_degree);
			auxiliary_points2.resize(v_count + 2 * m_v_degree);

            m_v_knots_vector = surf->get_v_knots();
            m_u_knots_vector = surf->get_u_knots();

			valid_auxiliary_u_knots_count = m_u_knots_vector.rows();
			valid_auxiliary_v_knots_count = m_v_knots_vector.rows();
			auxiliary_u_knots_vector.resize(valid_auxiliary_u_knots_count + 2 * m_u_degree);
			auxiliary_v_knots_vector.resize(valid_auxiliary_v_knots_count + 2 * m_v_degree);
			for (int index = 0; index < v_count + 2 * m_v_degree; ++index)
			{
				auxiliary_points[index].resize(point_size, u_count + 2 * m_u_degree);
				auxiliary_points2[index].resize(point_size, u_count + 2 * m_u_degree);
			}
			int degree = std::max(m_u_degree, m_v_degree);
			u_ders_basis_funs_array.resize(m_u_degree + 1, 3);
			v_ders_basis_funs_array.resize(m_v_degree + 1, 3);
			ndu.resize(degree + 1, degree + 1);
			ndu_trans.resize(degree + 1, degree + 1);
			iterArray.resize(degree + 2);

			left.resize(degree);
			right.resize(degree);
			coeff_denominator.resize(degree, 1);
			first_coeff.resize(degree + 1, 1);
			second_coeff.resize(degree + 1, 1);
			gammas.resize(degree + 1);
			gammas[0] = 1;
			for (int index = 1; index <= degree; ++index)
			{
				gammas[index] = index * gammas[index - 1];
			}
			arrays[0].resize(degree + 1, 1);
			arrays[1].resize(degree + 1, 1);

			temp.resize(point_size, m_v_degree + 1);
			RW.resize(point_size, degree + 1);
			RW2.resize(degree + 1);
			for (int index = 0; index <= degree; ++index)
			{
				RW2[index].resize(point_size, u_count + 2 * m_u_degree);
			}
			u_basis_funs_array.resize(m_u_degree + 1);
			v_basis_funs_array.resize(m_v_degree + 1);
		}

		~surf_compute() {};

		template<typename T, bool flag = ENUM_LIMITDIRECTION::RIGHT>
		ENUM_NURBS find_span(T u, int degree, Eigen::Vector<T, Eigen::Dynamic> const& knots_vector, int low, int high, int& index)
		{
			static_assert(flag == ENUM_LIMITDIRECTION::RIGHT || flag == ENUM_LIMITDIRECTION::LEFT, "find_span : flag have to equal RIGHT or LEFT");
			int points_count = high - low + degree;
			if (u > knots_vector[high] || u < knots_vector[low])
				return ENUM_NURBS::NURBS_PARAM_IS_OUT_OF_DOMAIN;
			int mid = (low + high) / 2;
			if constexpr (flag == ENUM_LIMITDIRECTION::RIGHT)
			{
				if (u == knots_vector[points_count])
				{
					index = points_count - 1;
					return ENUM_NURBS::NURBS_SUCCESS;
				}
				while (u < knots_vector[mid] || u >= knots_vector[mid + 1])
				{
					if (u < knots_vector[mid])
						high = mid;
					else
						low = mid;
					mid = (low + high) / 2;
				}
			}
			else
			{
				if (u == knots_vector[low])
				{
					index = low;
					return ENUM_NURBS::NURBS_SUCCESS;
				}

				while (u <= knots_vector[mid] || u > knots_vector[mid + 1])
				{
					if (u <= knots_vector[mid])
						high = mid;
					else
						low = mid;
					mid = (low + high) / 2;
				}
			}

			index = mid;
			return ENUM_NURBS::NURBS_SUCCESS;
		}
		ENUM_NURBS refine_knots_vector(const Eigen::VectorX<T>& insert_knots, ENUM_DIRECTION direction)
		{
			if (insert_knots.size() == 0)
				return ENUM_NURBS::NURBS_SUCCESS;

			if (direction == ENUM_DIRECTION::U_DIRECTION)
			{
				int cols = valid_auxiliary_u_knots_count - m_u_degree - 1;
				int rows = valid_auxiliary_v_knots_count - m_v_degree - 1;

				int start_index = -1, end_index = -1;
				if (find_span<T>(insert_knots[0], m_u_degree, auxiliary_u_knots_vector, 0, valid_auxiliary_u_knots_count - m_u_degree - 1, start_index) != ENUM_NURBS::NURBS_SUCCESS)
					return ENUM_NURBS::NURBS_ERROR;
				int insert_knots_size = insert_knots.size();
				if (find_span<T>(insert_knots[insert_knots_size - 1], m_u_degree, auxiliary_u_knots_vector, 0, valid_auxiliary_u_knots_count - m_u_degree - 1, end_index) != ENUM_NURBS::NURBS_SUCCESS)
					return ENUM_NURBS::NURBS_ERROR;
				end_index += 1;


				auxiliary_u_knots_vector.block(end_index + m_u_degree + insert_knots_size, 0, cols - end_index + 1, 1) = auxiliary_u_knots_vector.block(end_index + m_u_degree, 0, cols - end_index + 1, 1).eval();

				int i = end_index + m_u_degree - 1, k = end_index + insert_knots_size + m_u_degree - 1;
				for (int j = insert_knots_size - 1; j >= 0; --j)
				{
					while (insert_knots[j] <= auxiliary_u_knots_vector[i] && i > start_index)
					{
						auxiliary_u_knots_vector[k] = m_u_knots_vector[i];
						--k; --i;
					}
					auxiliary_u_knots_vector[k] = insert_knots[j];
					--k;
				}

				for (int row_index = 0; row_index < rows; ++row_index)
				{
					for (int index = cols - end_index; index >= 0; --index)
					{
						auxiliary_points[row_index].col(index + end_index + insert_knots_size - 1) = auxiliary_points[row_index].col(index + end_index - 1);
						// new_control_points.block(0, index + end_index + insert_knots_size - 1, point_size, 1) = m_control_points[row_index].col(index + end_index - 1);
					}

					i = end_index + m_u_degree - 1;
					k = end_index + insert_knots_size + m_u_degree - 1;
					for (int j = insert_knots_size - 1; j >= 0; --j)
					{
						while (insert_knots[j] <= m_u_knots_vector[i] && i > start_index)
						{
							auxiliary_points[row_index].col(k - m_u_degree - 1) = auxiliary_points[row_index].col(i - m_u_degree - 1);
							// new_control_points.block(0, k - m_u_degree - 1, point_size, 1) = m_control_points[row_index].col(i - m_u_degree - 1);
							--k; --i;
						}
						// new_control_points.block(0, k - m_u_degree - 1, point_size, 1) = new_control_points.block(0, k - m_u_degree, point_size, 1);
						auxiliary_points[row_index].col(k - m_u_degree - 1) = auxiliary_points[row_index].col(k - m_u_degree);
						for (int l = 1; l <= m_u_degree; ++l)
						{
							int ind = k - m_u_degree + l;
							T alpha = auxiliary_u_knots_vector[k + l] - insert_knots[j];
							if (alpha == 0.0)
								// new_control_points[row_index].col(ind - 1) = new_control_points[row_index].col(ind);
								auxiliary_points[row_index].col(ind - 1) = auxiliary_points[row_index].col(ind);
								// temp_control_points<T, point_size>::new_control_points.block(0, ind - 1, point_size, 1) = temp_control_points<T, point_size>::new_control_points.block(0, ind, point_size, 1);
							else
							{
								alpha /= (auxiliary_u_knots_vector[k + l] - m_u_knots_vector[i - m_u_degree + l]);
								// new_control_points.block(0, ind - 1, point_size, 1) *= alpha;
								// new_control_points.block(0, ind - 1, point_size, 1) += (1.0 - alpha) * new_control_points.block(0, ind, point_size, 1);
								auxiliary_points[row_index].col(ind - 1) *= alpha;
								auxiliary_points[row_index].col(ind - 1) += (1.0 - alpha) * auxiliary_points[row_index].col(ind);
							}
						}
						--k;
					}
				}
				valid_auxiliary_u_knots_count += insert_knots.rows();
			}
			else
			{
				int cols = valid_auxiliary_u_knots_count - m_u_degree - 1;
				int rows = valid_auxiliary_v_knots_count - m_v_degree - 1;

				int start_index = -1, end_index = -1;
				if (find_span<T>(insert_knots[0], m_v_degree, auxiliary_v_knots_vector, 0, valid_auxiliary_v_knots_count - m_v_degree - 1, start_index) != ENUM_NURBS::NURBS_SUCCESS)
					return ENUM_NURBS::NURBS_ERROR;
				int insert_knots_size = insert_knots.size();
				if (find_span<T>(insert_knots[insert_knots_size - 1], m_v_degree, auxiliary_v_knots_vector, 0, valid_auxiliary_v_knots_count - m_v_degree - 1, end_index) != ENUM_NURBS::NURBS_SUCCESS)
					return ENUM_NURBS::NURBS_ERROR;
				end_index += 1;

				for (int index = rows - end_index; index >= 0; --index)
				{
					auxiliary_points[index + end_index + insert_knots_size - 1].block(0, 0, point_size, cols)= auxiliary_points[index + end_index - 1].block(0, 0, point_size, cols);
				}

				auxiliary_v_knots_vector.block(end_index + m_v_degree + insert_knots_size, 0, rows - end_index + 1, 1) = auxiliary_v_knots_vector.block(end_index + m_v_degree, 0, rows - end_index + 1, 1).eval();

				int i = end_index + m_v_degree - 1, k = end_index + insert_knots_size + m_v_degree - 1;
				for (int j = insert_knots_size - 1; j >= 0; --j)
				{
					while (insert_knots[j] <= m_v_knots_vector[i] && i > start_index)
					{
						auxiliary_v_knots_vector[k] = m_v_knots_vector[i];
						--k; --i;
					}
					auxiliary_v_knots_vector[k] = insert_knots[j];
					--k;
				}

				i = end_index + m_v_degree - 1;
				k = end_index + insert_knots_size + m_v_degree - 1;
				for (int j = insert_knots_size - 1; j >= 0; --j)
				{
					while (insert_knots[j] <= m_v_knots_vector[i] && i > start_index)
					{
						auxiliary_points[k - m_v_degree - 1].block(0, 0, point_size, cols) = auxiliary_points[i - m_v_degree - 1].block(0, 0, point_size, cols);
						--k; --i;
					}
					// auxiliary_v_knots_vector[k - m_v_degree - 1] = auxiliary_v_knots_vector[k - m_v_degree];
					auxiliary_points[k - m_v_degree - 1].block(0, 0, point_size, cols) = auxiliary_points[k - m_v_degree].block(0, 0, point_size, cols);
					for (int l = 1; l <= m_v_degree; ++l)
					{
						int ind = k - m_v_degree + l;
						T alpha = auxiliary_v_knots_vector[k + l] - insert_knots[j];
						// if (std::abs(alpha) < KNOTS_VECTOR_EPS)
						if (alpha == 0.0)
						{
						    auxiliary_points[ind - 1].block(0, 0, point_size, cols) = auxiliary_points[ind].block(0, 0, point_size, cols);
						    
							// auxiliary_v_knots_vector[ind - 1] = auxiliary_v_knots_vector[ind];
						}
						else
						{
							alpha /= (auxiliary_v_knots_vector[k + l] - auxiliary_v_knots_vector[i - m_v_degree + l]);
							auxiliary_points[ind - 1].block(0, 0, point_size, cols) *= alpha;
							auxiliary_points[ind - 1].block(0, 0, point_size, cols) += (1.0 - alpha) * auxiliary_points[ind].block(0, 0, point_size, cols);
						}
					}
					--k;
				}
				valid_auxiliary_v_knots_count += insert_knots.rows();
			}

			return ENUM_NURBS::NURBS_SUCCESS;

		}
	
		ENUM_NURBS insert_knots(Box<T, 2>& sub_box, Box<int, 2>& rs)
		{
			int rows = m_v_knots_vector.rows() - m_v_degree - 1;
			const auto& old_control_points = m_surface->get_control_points_ref();
			int start_span = -1, end_span = -1;
			if (find_span<T>(sub_box.Min[0], m_u_degree, m_u_knots_vector, m_u_degree, m_u_knots_vector.rows() - m_u_degree - 1, start_span) != ENUM_NURBS::NURBS_SUCCESS)
				return ENUM_NURBS::NURBS_ERROR;

			if (find_span<T>(sub_box.Max[0], m_u_degree, m_u_knots_vector, m_u_degree, m_u_knots_vector.rows() - m_u_degree - 1, end_span) != ENUM_NURBS::NURBS_SUCCESS)
				return ENUM_NURBS::NURBS_ERROR;
			for (int index = 0; index < m_u_degree + 1; ++index)
			{
				auxiliary_u_knots_vector[index] = sub_box.Min[0];
			}
			for (int index = start_span + 1; index <= end_span + 1 + m_u_degree; ++index)
			{
				auxiliary_u_knots_vector[index - start_span - 1 + m_u_degree + 1] = m_u_knots_vector[index];
			}
			// end_span += r1;
			//插入节点r次
			for (int index = 0; index < rows; ++index)
			{
				int repeat_count = m_u_degree - rs.Min[0];
				RW = old_control_points[index].block(0, start_span - m_u_degree, point_size, m_u_degree - repeat_count + 1);
				for (int j = 1; j <= rs.Min[0]; ++j)
				{
					int L = start_span - m_u_degree + j;
					int count = m_u_degree - j - repeat_count;
					for (int i = 0; i <= count; ++i)
					{
						T alpha = (sub_box.Min[0] - m_u_knots_vector[i + L]) / (m_u_knots_vector[i + start_span + 1] - m_u_knots_vector[i + L]);
						RW.col(i) = alpha * RW.col(i + 1) + (1.0 - alpha) * RW.col(i);
					}
					auxiliary_points[index].col(rs.Min[0] - j) = RW.col(m_u_degree - j - repeat_count);
				}

				int new_span = end_span - start_span + rs.Min[0];
				new_span = std::max(m_u_degree, new_span);
				valid_auxiliary_u_knots_count = new_span + m_u_degree + 2;
				auxiliary_points[index].block(0, rs.Min[0], point_size, repeat_count + 1) = old_control_points[index].block(0, start_span - repeat_count, point_size, repeat_count + 1);
				// for (int j = new_span; index < new_span + repeat_count + 1; ++index)
				// {
				// 	auxiliary_points2[index].block(0, 0, point_size, cols) = auxiliary_points[start_span + (index - new_span)].block(0, 0, point_size, cols);
				// }
				repeat_count = m_u_degree - rs.Max[0];
				RW = auxiliary_points[index].block(0, new_span - m_u_degree, point_size, m_u_degree - repeat_count + 1);
				for (int j = 1; j <= rs.Max[0]; ++j)
				{
					int L = new_span - m_u_degree + j;
					int count = m_u_degree - j - repeat_count;
					for (int i = 0; i <= count; ++i)
					{
						T alpha = (sub_box.Max[0] - auxiliary_u_knots_vector[i + L]) / (auxiliary_u_knots_vector[i + new_span + 1] - auxiliary_u_knots_vector[i + L]);
						RW.col(i) = alpha * RW.col(i + 1) + (1.0 - alpha) * RW.col(i);
					}
					auxiliary_points[index].col(L) = RW.col(0);
					// auxiliary_points[index].col(new_span + rs.Max[0] - j - repeat_count) = RW.col(m_u_degree - j - repeat_count);
				}
			}
			for (int index = 0; index < m_u_degree + 1; ++index)
			{
				auxiliary_u_knots_vector[index + end_span - start_span + m_u_degree + 1] = sub_box.Max[0];
			}

			// valid_auxiliary_u_knots_count = end_span - start_span + rs.Min[0] + rs.Max[0] + 2;
			// valid_auxiliary_u_knots_count = 2 * m_u_degree + 2;
			// std::cout << "begin1" << std::endl;
			// for (int index = 0; index < auxiliary_points.size(); ++index)
			// {
			// 	std::cout << auxiliary_points[index] << std::endl;
			// 	std::cout << "****************************" << std::endl;
			// }
			// std::cout << "end1" << std::endl;

			int cols = valid_auxiliary_u_knots_count - m_u_degree - 1;

			start_span = -1; 
			end_span = -1;
			if (find_span<T>(sub_box.Min[1], m_v_degree, m_v_knots_vector, m_v_degree, m_v_knots_vector.rows() - m_v_degree - 1, start_span) != ENUM_NURBS::NURBS_SUCCESS)
				return ENUM_NURBS::NURBS_ERROR;
			if (find_span<T>(sub_box.Max[1], m_v_degree, m_v_knots_vector, m_v_degree, m_v_knots_vector.rows() - m_v_degree - 1, end_span) != ENUM_NURBS::NURBS_SUCCESS)
				return ENUM_NURBS::NURBS_ERROR;
			// end_index += 1;

			for (int index = 0; index < m_v_degree + 1; ++index)
			{
				auxiliary_v_knots_vector[index] = sub_box.Min[1];
			}
			for (int index = start_span + 1; index <= end_span + 1 + m_v_degree; ++index)
			{
				auxiliary_v_knots_vector[index - start_span - 1 + m_v_degree + 1] = m_v_knots_vector[index];
				// u_temp_knots.push_back(m_u_knots_vector[index]);
			}


			//插入节点r次
			int repeat_count = m_v_degree - rs.Min[1];
			for (int index = 0; index < m_v_degree - repeat_count + 1; ++index)
			{
				RW2[index].block(0, 0, point_size, cols) = auxiliary_points[index + start_span - m_v_degree].block(0, 0, point_size, cols);
			}
			// RW = old_control_points[index].block(0, start_span - m_u_degree, point_size, m_v_degree - repeat_count + 1);
			for (int j = 1; j <= rs.Min[1]; ++j)
			{
				int L = start_span - m_v_degree + j;
				int count = m_v_degree - j - repeat_count;
				for (int i = 0; i <= count; ++i)
				{
					T alpha = (sub_box.Min[1] - m_v_knots_vector[i + L]) / (m_v_knots_vector[i + start_span + 1] - m_v_knots_vector[i + L]);
					RW2[i].block(0, 0, point_size, cols) = alpha * RW2[i + 1].block(0, 0, point_size, cols) + (1.0 - alpha) * RW2[i].block(0, 0, point_size, cols);
				}
			//	old_control_points[index].block(0, start_span, point_size, end_span - start_span);
				auxiliary_points2[rs.Min[1] - j].block(0, 0, point_size, cols) = RW2[m_v_degree - j - repeat_count].block(0, 0, point_size, cols);
			}
			int new_span = end_span - start_span + rs.Min[1];
			valid_auxiliary_v_knots_count = new_span + m_v_degree + 2;
			new_span = std::max(m_v_degree, new_span);
			for (int index = 0; index < repeat_count + 1; ++index)
			{
				auxiliary_points2[rs.Min[1] + index].block(0, 0, point_size, cols) = auxiliary_points[start_span + index - repeat_count].block(0, 0, point_size, cols);
			}
			repeat_count = m_v_degree - rs.Max[1];
			for (int index = 0; index < m_v_degree - repeat_count + 1; ++index)
			{
				RW2[index].block(0, 0, point_size, cols) = auxiliary_points2[index + new_span - m_v_degree].block(0, 0, point_size, cols);
			}
			for (int j = 1; j <= rs.Max[1]; ++j)
			{
				int L = new_span - m_v_degree + j;
				int count = m_v_degree - j - repeat_count;
				for (int i = 0; i <= count; ++i)
				{
					T alpha = (sub_box.Max[1] - auxiliary_v_knots_vector[i + L]) / (auxiliary_v_knots_vector[i + new_span + 1] - auxiliary_v_knots_vector[i + L]);
					RW2[i].block(0, 0, point_size, cols) = alpha * RW2[i + 1].block(0, 0, point_size, cols) + (1.0 - alpha) * RW2[i].block(0, 0, point_size, cols);
				}
				auxiliary_points2[L].block(0, 0, point_size, cols) = RW2[0].block(0, 0, point_size, cols);
			}
			for (int index = 0; index < m_v_degree + 1; ++index)
			{
				auxiliary_v_knots_vector[index + end_span - start_span + m_v_degree + 1] = sub_box.Max[1];
			}
			return ENUM_NURBS::NURBS_SUCCESS;

		}

		ENUM_NURBS split_at_param(surface_type* surf ,const Eigen::Vector2<T>& param, std::vector<surface_type*>& sub_nurbs)
		{
			std::array<int, 2> rs{ m_u_degree, m_v_degree };
			const auto& old_control_points = surf->get_control_points_ref();
			const Eigen::VectorX<T>& u_knots_vector = surf->m_u_knots_vector;
			const Eigen::VectorX<T>& v_knots_vector = surf->m_v_knots_vector;
			int rows = v_knots_vector.rows() - m_v_degree - 1;
			int start_span = -1;
			if (find_span<T>(param[0], m_u_degree, u_knots_vector, m_u_degree, u_knots_vector.rows() - m_u_degree - 1, start_span) != ENUM_NURBS::NURBS_SUCCESS)
				return ENUM_NURBS::NURBS_ERROR;

			if (param[0] == u_knots_vector[start_span])
			{
				rs[0] = 1;
				int index = start_span - 1;
				while (index >= 0)
				{
					if (param[0] == u_knots_vector[index])
					{
						--index;
						--rs[0];
					}
					else
					{
						break;
					}
				}
			}
			else if (param[0] == u_knots_vector[start_span + 1])
			{
				rs[0] = 1;
				int index = start_span + 2;
				while (index < u_knots_vector.rows())
				{
					if (param[0] == u_knots_vector[index])
					{
						++index;
						--rs[0];
					}
					else
					{
						break;
					}
				}
			}
			rs[0] = std::max(0, rs[0]);

			//插入节点r次
			int current_cols_count = u_knots_vector.rows() - m_u_degree - 1;
			int repeat_count = m_u_degree - rs[0];
			for (int index = 0; index < rows; ++index)
			{
				auxiliary_points[index].middleCols(0, start_span - m_u_degree + 1) = old_control_points[index].middleCols(0, start_span - m_u_degree + 1);
				auxiliary_points[index].middleCols(start_span - repeat_count + rs[0] + 1, current_cols_count - start_span + repeat_count) = old_control_points[index].middleCols(start_span - repeat_count, current_cols_count - start_span + repeat_count);
				RW = old_control_points[index].block(0, start_span - m_u_degree, point_size, m_u_degree - repeat_count + 1);
				for (int j = 1; j <= rs[0]; ++j)
				{
					int L = start_span - m_u_degree + j;
					int count = m_u_degree - j - repeat_count;
					for (int i = 0; i <= count; ++i)
					{
						T alpha = (param[0] - u_knots_vector[i + L]) / (u_knots_vector[i + start_span + 1] - u_knots_vector[i + L]);
						RW.col(i) = alpha * RW.col(i + 1) + (1.0 - alpha) * RW.col(i);
					}
					auxiliary_points[index].col(L) = RW.col(0);
					auxiliary_points[index].col(rs[0] - j + start_span - repeat_count + 1) = RW.col(m_u_degree - j - repeat_count);
				}
			}

			int cols = u_knots_vector.rows() - m_u_degree + rs[0];
			int u_start_span = start_span;
			int u_repeat_count = repeat_count;

			start_span = -1; 
			if (find_span<T>(param[1], m_v_degree, v_knots_vector, m_v_degree, v_knots_vector.rows() - m_v_degree - 1, start_span) != ENUM_NURBS::NURBS_SUCCESS)
				return ENUM_NURBS::NURBS_ERROR;

			if (param[1] == v_knots_vector[start_span])
			{
				rs[1] = 1;
				int index = start_span - 1;
				while (index >= 0)
				{
					if (param[1] == v_knots_vector[index])
					{
						--index;
						--rs[1];
					}
					else
					{
						break;
					}
				}
			}
			else if (param[1] == v_knots_vector[start_span + 1])
			{
				rs[1] = 1;
				int index = start_span + 2;
				while (index < v_knots_vector.rows())
				{
					if (param[1] == v_knots_vector[index])
					{
						++index;
						--rs[1];
					}
					else
					{
						break;
					}
				}
			}
			rs[1] = std::max(0, rs[1]);
			//插入节点r次
			repeat_count = m_v_degree - rs[1];
			for (int index = 0; index < m_v_degree - repeat_count + 1; ++index)
			{
				RW2[index].block(0, 0, point_size, cols) = auxiliary_points[index + start_span - m_v_degree].block(0, 0, point_size, cols);
			}

			// printNurbsControlPoints(auxiliary_points);
			sub_nurbs.resize(4);
			for (int index = 0; index < 2; ++index)
			{
				sub_nurbs[index] = new surface_type();
				sub_nurbs[index]->m_control_points.resize(start_span + 1 - repeat_count);
				sub_nurbs[index]->set_uv_degree(m_u_degree, m_v_degree);
			}
			for (int i = 0; i < start_span + 1 - m_v_degree; ++i)
			{
				sub_nurbs[0]->m_control_points[i] = std::move(auxiliary_points[i].middleCols(0, u_start_span + 1 - u_repeat_count));
				sub_nurbs[1]->m_control_points[i] = std::move(auxiliary_points[i].middleCols(u_start_span + 1 - u_repeat_count, u_knots_vector.rows() - u_start_span - 1));
			}
			// printNurbsControlPoints(sub_nurbs[1]->m_control_points);


			for (int index = 2; index < 4; ++index)
			{
				sub_nurbs[index] = new surface_type();
				sub_nurbs[index]->m_control_points.resize(v_knots_vector.rows() - start_span - 1);
				sub_nurbs[index]->set_uv_degree(m_u_degree, m_v_degree);
			}
			for (int i = 0; i < rows - start_span + repeat_count; ++i)
			{
				sub_nurbs[2]->m_control_points[i + (m_v_degree - repeat_count)] = std::move(auxiliary_points[start_span - repeat_count + i].middleCols(0, u_start_span + 1 - u_repeat_count));
				sub_nurbs[3]->m_control_points[i + (m_v_degree - repeat_count)] = std::move(auxiliary_points[start_span - repeat_count + i].middleCols(u_start_span + 1 - u_repeat_count, u_knots_vector.rows() - u_start_span - 1));
			}

			sub_nurbs[0]->m_u_knots_vector = u_knots_vector.head(u_start_span + 1 + (m_u_degree + 1) - u_repeat_count);
			sub_nurbs[0]->m_u_knots_vector.tail(m_u_degree + 1).setConstant(param[0]);

			sub_nurbs[0]->m_v_knots_vector = v_knots_vector.head(start_span + 1 + (m_v_degree + 1) - repeat_count);
			sub_nurbs[0]->m_v_knots_vector.tail(m_v_degree + 1).setConstant(param[1]);

			sub_nurbs[1]->m_u_knots_vector = u_knots_vector.tail(u_knots_vector.rows() - u_start_span + m_u_degree);
			sub_nurbs[1]->m_u_knots_vector.head(m_u_degree + 1).setConstant(param[0]);

			sub_nurbs[1]->m_v_knots_vector = v_knots_vector.head(start_span + 1 + (m_v_degree + 1) - repeat_count);
			sub_nurbs[1]->m_v_knots_vector.tail(m_v_degree + 1).setConstant(param[1]);
			
			sub_nurbs[2]->m_u_knots_vector = u_knots_vector.head(u_start_span + 1 + (m_u_degree + 1) - u_repeat_count);
			sub_nurbs[2]->m_u_knots_vector.tail(m_u_degree + 1).setConstant(param[0]);

			sub_nurbs[2]->m_v_knots_vector = v_knots_vector.tail(v_knots_vector.rows() - start_span + m_v_degree);
			sub_nurbs[2]->m_v_knots_vector.head(m_v_degree + 1).setConstant(param[1]);
			
			sub_nurbs[3]->m_u_knots_vector = u_knots_vector.tail(u_knots_vector.rows() - u_start_span + m_u_degree);
			sub_nurbs[3]->m_u_knots_vector.head(m_u_degree + 1).setConstant(param[0]);

			sub_nurbs[3]->m_v_knots_vector = v_knots_vector.tail(v_knots_vector.rows() - start_span + m_v_degree);
			sub_nurbs[3]->m_v_knots_vector.head(m_v_degree + 1).setConstant(param[1]);
			
			for (int j = 1; j <= rs[1]; ++j)
			{
				int L = start_span - m_v_degree + j;
				int count = m_v_degree - j - repeat_count;
				for (int i = 0; i <= count; ++i)
				{
					T alpha = (param[1] - v_knots_vector[i + L]) / (v_knots_vector[i + start_span + 1] - v_knots_vector[i + L]);
					RW2[i].block(0, 0, point_size, cols) = alpha * RW2[i + 1].block(0, 0, point_size, cols) + (1.0 - alpha) * RW2[i].block(0, 0, point_size, cols);
				}

				sub_nurbs[0]->m_control_points[L] = RW2[0].middleCols(0, u_start_span + 1 - u_repeat_count);
				sub_nurbs[1]->m_control_points[L] = RW2[0].middleCols(u_start_span + 1 - u_repeat_count, u_knots_vector.rows() - u_start_span - 1);
				sub_nurbs[2]->m_control_points[rs[1] - j] = RW2[m_v_degree - j - repeat_count].middleCols(0, u_start_span + 1 - u_repeat_count);
				sub_nurbs[3]->m_control_points[rs[1] - j] = RW2[m_v_degree - j - repeat_count].middleCols(u_start_span + 1 - u_repeat_count, u_knots_vector.rows() - u_start_span - 1);
			}
			// for (int i = 0; i < 4; ++i)
			// {
			// 	printNurbsControlPoints(sub_nurbs[i]->m_control_points);
			// }

			return ENUM_NURBS::NURBS_SUCCESS;

		}
		void get_2_ders_sub_box2(Box<T, 2>& uv_box, Box<T, dim>& u_tangent_box, Box<T, dim>& v_tangent_box)
		{
			static_assert(is_rational == false, "is_ratio != false");

			Box<int, 2> rs;

			int u_begin_index;
			int u_begin_mul = konts_multiple<T>(uv_box.Min[0], m_u_knots_vector, m_u_degree, u_begin_index);
			rs.Min[0] = std::max(0, m_u_degree - u_begin_mul);
			int u_end_index;
			int u_end_mul = konts_multiple<T>(uv_box.Max[0], m_u_knots_vector, m_u_degree, u_end_index);
			rs.Max[0] = std::max(m_u_degree - u_end_mul, 0);

			int v_begin_index;
			int v_begin_mul = konts_multiple<T>(uv_box.Min[1], m_v_knots_vector, m_v_degree, v_begin_index);
			rs.Min[1] = std::max(m_v_degree - v_begin_mul, 0);
			int v_end_index;
			int v_end_mul = konts_multiple<T>(uv_box.Max[1], m_v_knots_vector, m_v_degree, v_end_index);
			rs.Max[1] = std::max(m_v_degree - v_end_mul, 0);

			insert_knots(uv_box, rs);

			Eigen::Vector<T, dim> point = (m_u_degree / (auxiliary_u_knots_vector[1 + m_u_degree] - auxiliary_u_knots_vector[1])) * Eigen::Vector<T, dim>((auxiliary_points2[0].col(1) - auxiliary_points2[0].col(0)));
			u_tangent_box.Min = u_tangent_box.Max = point;
			for (Eigen::Index u_index = 0; u_index < valid_auxiliary_u_knots_count - m_u_degree - 2; ++u_index)
			{
				T coeff = m_u_degree / (auxiliary_u_knots_vector[u_index + 1 + m_u_degree] - auxiliary_u_knots_vector[u_index + 1]);
				for (Eigen::Index v_index = 0; v_index < valid_auxiliary_v_knots_count - m_v_degree - 1; ++v_index)
				{
					point = coeff * Eigen::Vector<T, dim>((auxiliary_points2[v_index].col(1 + u_index) - auxiliary_points2[v_index].col(u_index)));
					u_tangent_box.enlarge(point);
				}
			}

			point = (m_v_degree / (auxiliary_v_knots_vector[1 + m_v_degree] - auxiliary_v_knots_vector[1])) * Eigen::Vector<T, dim>((auxiliary_points2[1].col(0) - auxiliary_points2[0].col(0)));
			v_tangent_box.Min = v_tangent_box.Max = point;
			for (Eigen::Index v_index = 0; v_index < valid_auxiliary_v_knots_count - m_v_degree - 2; ++v_index)
			{
				T coeff = (m_v_degree / (auxiliary_v_knots_vector[v_index + 1 + m_v_degree] - auxiliary_v_knots_vector[v_index + 1]));
				for (Eigen::Index u_index = 0; u_index < valid_auxiliary_u_knots_count - m_u_degree - 1; ++u_index)
				{
					point = coeff * Eigen::Vector<T, dim>((auxiliary_points2[v_index + 1].col(u_index) - auxiliary_points2[v_index].col(u_index)));
					v_tangent_box.enlarge(point);
				}
			}

			return;
		}

		void get_2_ders_sub_box(Box<T, 2>& uv_box, Box<T, dim>& u_tangent_box, Box<T, dim>& v_tangent_box)
        {
            static_assert(is_rational == false, "is_ratio != false");
			
			const auto& control_points = m_surface->get_control_points_ref();
			int v_count = control_points.size();
			int u_count = control_points[0].cols();

			for (int index = 0; index < v_count; ++index)
			{
				auxiliary_points[index].block(0, 0, point_size, u_count) = control_points[index];
			}
			auxiliary_u_knots_vector.block(0, 0, m_u_knots_vector.rows(), 1) = m_u_knots_vector;
			auxiliary_v_knots_vector.block(0, 0, m_v_knots_vector.rows(), 1) = m_v_knots_vector;
			valid_auxiliary_u_knots_count = m_u_knots_vector.rows();
			valid_auxiliary_v_knots_count = m_v_knots_vector.rows();

            int u_begin_index;
            int u_begin_mul = konts_multiple<T>(uv_box.Min[0], m_u_knots_vector, m_u_degree, u_begin_index);
            int u_begin_insert_num = std::max(0, m_u_degree - u_begin_mul);
            int u_end_index;
            int u_end_mul = konts_multiple<T>(uv_box.Max[0], m_u_knots_vector, m_u_degree, u_end_index);
            int u_end_insert_num = std::max(m_u_degree - u_end_mul, 0);
            Eigen::VectorX<T> insert_knots(u_begin_insert_num + u_end_insert_num);
            insert_knots.block(0, 0, u_begin_insert_num, 1).setConstant(uv_box.Min[0]);
            insert_knots.block(u_begin_insert_num, 0, u_end_insert_num, 1).setConstant(uv_box.Max[0]);
            refine_knots_vector(insert_knots, ENUM_DIRECTION::U_DIRECTION);
			std::cout << "begin1" << std::endl;
			for (int index = 0; index < auxiliary_points.size(); ++index)
			{
				std::cout << auxiliary_points[index] << std::endl;
				std::cout << "****************************" << std::endl;
			}
			std::cout << "end1" << std::endl;

            int v_begin_index;
            int v_begin_mul = konts_multiple<T>(uv_box.Min[1], m_v_knots_vector, m_v_degree, v_begin_index);
            int v_begin_insert_num = std::max(m_v_degree - v_begin_mul, 0);
            int v_end_index;
            int v_end_mul = konts_multiple<T>(uv_box.Max[1], m_v_knots_vector, m_v_degree, v_end_index);
            int v_end_insert_num = std::max(m_v_degree - v_end_mul, 0);
            insert_knots.resize(v_begin_insert_num + v_end_insert_num);
            insert_knots.block(0, 0, v_begin_insert_num, 1).setConstant(uv_box.Min[1]);
            insert_knots.block(v_begin_insert_num, 0, v_end_insert_num, 1).setConstant(uv_box.Max[1]);
            refine_knots_vector(insert_knots, ENUM_DIRECTION::V_DIRECTION);
			std::cout << "begin2" << std::endl;
			for (int index = 0; index < auxiliary_points.size(); ++index)
			{
				std::cout << auxiliary_points[index] << std::endl;
				std::cout << "****************************" << std::endl;
			}
			std::cout << "end2" << std::endl;
            int u_new_knots_count = 2 * m_u_degree + 2 + u_end_index - u_begin_index - u_begin_mul;
            int v_new_knots_count = 2 * m_v_degree + 2 + v_end_index - v_begin_index - v_begin_mul;
     
			int v_start = std::max(v_begin_index - 1, 0);
            int u_start = std::max(u_begin_index - 1, 0);



            std::size_t row_points_count = v_new_knots_count - m_v_degree - 1;
            Eigen::Index col_points_count = u_new_knots_count - m_u_degree - 1;

            Eigen::Vector<T, dim> point = (m_u_degree / (auxiliary_u_knots_vector[1 + m_u_degree + u_start] - auxiliary_u_knots_vector[1 + u_start])) * Eigen::Vector<T, dim>((auxiliary_points[v_start].col(u_start + 1) - auxiliary_points[v_start].col(u_start)));
            u_tangent_box.Min = u_tangent_box.Max = point;
            for (Eigen::Index u_index = 0; u_index < col_points_count - 1; ++u_index)
            {
				T coeff = m_u_degree / (auxiliary_u_knots_vector[u_index + 1 + m_u_degree + u_start] - auxiliary_u_knots_vector[u_index + 1 + u_start]);
				for (Eigen::Index v_index = 0; v_index < row_points_count; ++v_index)
                {
					point = coeff * Eigen::Vector<T, dim>((auxiliary_points[v_index + v_start].col(u_start + 1 + u_index) - auxiliary_points[v_index + v_start].col(u_index + u_start)));
                    u_tangent_box.enlarge(point);
                }
            }

            point = (m_v_degree / (auxiliary_v_knots_vector[1 + m_v_degree + v_start] - auxiliary_v_knots_vector[1 + v_start])) * Eigen::Vector<T, dim>((auxiliary_points[v_start + 1].col(u_start) - auxiliary_points[v_start].col(u_start)));
            v_tangent_box.Min = v_tangent_box.Max = point;
            for (Eigen::Index v_index = 0; v_index < row_points_count - 1; ++v_index)
            {
				T coeff = (m_v_degree / (auxiliary_v_knots_vector[v_index + 1 + m_v_degree + v_start] - auxiliary_v_knots_vector[v_index + 1 + v_start]));
				for (Eigen::Index u_index = 0; u_index < col_points_count; ++u_index)
                {
					point = coeff * Eigen::Vector<T, dim>((auxiliary_points[v_index + 1 + v_start].col(u_index + u_start) - auxiliary_points[v_index + v_start].col(u_index + u_start)));
                    v_tangent_box.enlarge(point);
                }
            }

			std::cout << "u_tangent: " << std::endl;
			std::cout << u_tangent_box.Min << std::endl;
			std::cout << u_tangent_box.Max << std::endl;
			std::cout << "uend" << std::endl;
			std::cout << "v_tangent: " << std::endl;
			std::cout << v_tangent_box.Min << std::endl;
			std::cout << v_tangent_box.Max << std::endl;
			std::cout << "vend" << std::endl;
            return;
        }
		
		template<typename T>
		ENUM_NURBS basis_functions(int i, T u, int degree, const Eigen::VectorX<T>& knots_vector, Eigen::Vector<T, Eigen::Dynamic>& result)
		{
			// 可以将下面left和right的计算改成多线程, 不过节点向量应该不大, 该不该影响应该不大
			//将u-u_j(j = i - p + 1, ... , i, 理论上应该算到u-u_(i - p), 但是其对应的基函数为0, 因此可以设置为0)
			// Eigen::Vector<T, Eigen::Dynamic> left(degree);
			//将u_j-u(j = i + 1, ... , i + p, 理论上应该算到u_(i + p + 1)-u, 但是其对应的基函数为0, 因此可以设置为0)
			// Eigen::Vector<T, Eigen::Dynamic> right(degree);
			for (int index = 0; index < degree; ++index)
			{
				left[index] = u - knots_vector[i - degree + 1 + index];
				right[index] = knots_vector[i + 1 + index] - u;
			}
			iterArray.setConstant(0.0);
			iterArray[degree] = 1.0;
			for (int iter_step = 1; iter_step <= degree; ++iter_step)
			{
				//u - u_j 列
				// Eigen::Array<T, Eigen::Dynamic, 1> first_coeff(iter_step + 1, 1);
				first_coeff[0] = 0.0;
				// Eigen::Array<T, Eigen::Dynamic, 1> second_coeff(iter_step + 1, 1);
				second_coeff[iter_step] = 0.0;

				const auto& first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
				const auto& second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

				const auto& coeff_denominator = first_coeff_numerator + second_coeff_numerator;

				first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator;
				second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator;
				first_coeff.block(0, 0, iter_step + 1, 1) *= iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
				second_coeff.block(0, 0, iter_step + 1, 1) *= iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();
				iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = first_coeff.block(0, 0, iter_step + 1, 1) + second_coeff.block(0, 0, iter_step + 1, 1);
			}
			result.block(0, 0, degree + 1, 1) = std::move(iterArray.block(0, 0, degree + 1, 1));
			return ENUM_NURBS::NURBS_SUCCESS;
		}
		
		
		template<typename T>
		ENUM_NURBS ders_basis_funs(int i, int n, int degree, T u, const Eigen::Vector<T, Eigen::Dynamic>& knots_vector, Eigen::MatrixX<T>& ders_basis_funs_array)
		{
			int new_n = std::min(n, degree);
			// result.resize(degree + 1, n + 1);
			ders_basis_funs_array.setConstant(0.0);
			// result.setConstant(0.0);
			ndu.setConstant(0.0);
			for (int index = 0; index < degree; ++index)
			{
				left[index] = u - knots_vector[i - degree + 1 + index];
				right[index] = knots_vector[i + 1 + index] - u;
			}
			// Eigen::Vector<T, Eigen::Dynamic> iterArray(degree + 2);
			iterArray.setConstant(0.0);
			iterArray[degree] = 1.0;
			ndu(0, 0) = 1.0;

			for (int iter_step = 1; iter_step <= degree; ++iter_step)
			{
				first_coeff[0] = 0.0;
				second_coeff[iter_step] = 0.0;

				const auto& first_coeff_numerator = left.block(degree - iter_step, 0, iter_step, 1).array();
				const auto& second_coeff_numerator = right.block(0, 0, iter_step, 1).array();

				coeff_denominator.block(0, 0, iter_step, 1) = first_coeff_numerator + second_coeff_numerator;

				first_coeff.block(1, 0, iter_step, 1) = first_coeff_numerator / coeff_denominator.block(0, 0, iter_step, 1);
				second_coeff.block(0, 0, iter_step, 1) = second_coeff_numerator / coeff_denominator.block(0, 0, iter_step, 1);
				first_coeff.block(0, 0, iter_step + 1, 1) *= iterArray.block(degree - iter_step, 0, iter_step + 1, 1).array();
				second_coeff.block(0, 0, iter_step + 1, 1) *= iterArray.block(degree - iter_step + 1, 0, iter_step + 1, 1).array();;

				iterArray.block(degree - iter_step, 0, iter_step + 1, 1) = first_coeff.block(0, 0, iter_step + 1, 1) + second_coeff.block(0, 0, iter_step + 1, 1);
				ndu.block(0, iter_step, iter_step + 1, 1) = iterArray.block(degree - iter_step, 0, iter_step + 1, 1);

				ndu.block(iter_step, 0, 1, iter_step) = coeff_denominator.block(0, 0, iter_step, 1).transpose();
			}

			ders_basis_funs_array.block(0, 0, degree + 1, 1) = ndu.block(0, degree, degree + 1, 1);

			ndu_trans = ndu.transpose();
			for (int r = i - degree; r <= i; ++r) //对基函数进行循环
			{
				arrays[0].setConstant(0.0);
				arrays[1].setConstant(0.0);
				int current_index = 0;
				int next_index = 1;

				arrays[0][0] = 1.0;
				for (int k = 1; k <= new_n; ++k)
				{
					auto& current_array = arrays[current_index];
					auto& next_array = arrays[next_index];
					next_array.setConstant(0.0);

					int left_num = r - (i - degree) - k;
					int left_index = std::max(0, -left_num);
					int right_num = r - (i - degree);
					int right_index = right_num > (degree - k) ? degree - k - left_num : right_num - left_num;
					int next_array_length = right_index - left_index + 1;

					int col_of_array = std::max(0, degree - i + r - k);
					int current_left_index = left_index;
					int current_right_index = right_index;
					if (left_index == 0)
					{
						next_array[0] = current_array[0] / ndu_trans(col_of_array, degree + 1 - k);
						current_left_index += 1;
					}
					if (right_index == k)
					{
						next_array[k] = -current_array[k - 1] / ndu_trans(next_array_length + col_of_array - 1, degree + 1 - k);
						current_right_index -= 1;
					}
					if (current_right_index >= current_left_index)
					{
						int arrayLength = current_right_index - current_left_index + 1;
						next_array.block(current_left_index, 0, arrayLength, 1) = (current_array.block(current_left_index, 0, arrayLength, 1) -
							current_array.block(current_left_index - 1, 0, arrayLength, 1)) / ndu_trans.block(col_of_array + current_left_index - left_index, degree + 1 - k, arrayLength, 1).array();
					}

					left_num = std::max(0, left_num);
					Eigen::Map<Eigen::VectorX<T>> temp_l(&ndu(left_num, degree - k), next_array_length);
					Eigen::Map<Eigen::VectorX<T>> temp_r(&next_array(left_index, 0), next_array_length);

					ders_basis_funs_array(r - (i - degree), k) = gammas[degree] / gammas[degree - k] *
						temp_l.dot(temp_r);
					std::swap(current_index, next_index);
				}
			}

			return ENUM_NURBS::NURBS_SUCCESS;
		}

		template<int n>
		ENUM_NURBS derivative_on_surface(T u, T v, Eigen::Matrix<Eigen::Vector<T, dim>, n + 1, n + 1>& result)
		{
			if constexpr (is_rational == false)
			{
				int du = std::min(n, m_u_degree);
				int dv = std::min(n, m_v_degree);
				// Eigen::Vector<T, point_size> zero_vector;
				// zero_vector.setConstant(0.0);
				// for (int k = m_u_degree + 1; k <= n; ++k)
				// {
				// 	int count = n - k;
				// 	for (int l = 0; l <= count; ++l)
				// 		result(k, l) = zero_vector;
				// }

				// for (int l = m_v_degree + 1; l <= n; ++l)
				// {
				// 	int count = n - l;
				// 	for (int k = 0; k <= count; ++k)
				// 		result(k, l) = zero_vector;
				// }

				int uspan = -1;
				int vspan = -1;       
				find_span<T>(u, m_u_degree, m_u_knots_vector, m_u_degree, m_u_knots_vector.rows() - m_u_degree - 1, uspan);
				find_span<T>(v, m_v_degree, m_v_knots_vector, m_v_degree, m_v_knots_vector.rows() - m_v_degree - 1, vspan);
				ders_basis_funs<T>(uspan, du, m_u_degree, u, m_u_knots_vector, u_ders_basis_funs_array);
				ders_basis_funs<T>(vspan, dv, m_v_degree, v, m_v_knots_vector, v_ders_basis_funs_array);
				const auto& control_points = m_surface->get_control_points_ref();
				for (int k = 0; k <= du; ++k)
				{
					temp.setConstant(0.0);
					for (int s = 0; s <= m_v_degree; ++s)
					{
						// temp.col(s) = control_points[vspan - m_v_degree + s].block(0, uspan - m_u_degree, point_size, m_u_degree + 1) * u_ders_basis_funs_array.col(k);
						temp.col(s) = control_points[vspan - m_v_degree + s].middleCols(uspan - m_u_degree, m_u_degree + 1) * u_ders_basis_funs_array.col(k);
					}
					int dd = std::min(n - k, dv);
					for (int l = 0; l <= dd; ++l)
					{
						result(k, l) = temp * v_ders_basis_funs_array.col(l);
					}
				}
				return ENUM_NURBS::NURBS_SUCCESS;
			}
			else
			{
				Eigen::Matrix<Eigen::Vector<T, point_size>, n + 1, n + 1> ration_result;
				int du = std::min(n, m_u_degree);
				int dv = std::min(n, m_v_degree);
				Eigen::Vector<T, point_size> zero_vector;
				zero_vector.setConstant(0.0);
				for (int k = m_u_degree + 1; k <= n; ++k)
				{
					int count = n - k;
					for (int l = 0; l <= count; ++l)
						ration_result(k, l) = zero_vector;
				}

				for (int l = m_v_degree + 1; l <= n; ++l)
				{
					int count = n - l;
					for (int k = 0; k <= count; ++k)
						ration_result(k, l) = zero_vector;
				}

				int uspan = -1;
				int vspan = -1;
				find_span<T>(u, m_u_degree, m_u_knots_vector, m_u_degree, m_u_knots_vector.rows() - m_u_degree - 1, uspan);
				find_span<T>(v, m_v_degree, m_v_knots_vector, m_u_degree, m_v_knots_vector.rows() - m_v_degree - 1, vspan);
				// find_span<T>(u, m_u_degree, m_u_knots_vector, uspan);
				// find_span<T>(v, m_v_degree, m_v_knots_vector, vspan);
				Eigen::MatrixX<T> nu, nv;
				ders_basis_funs<T>(uspan, du, m_u_degree, u, m_u_knots_vector, nu);
				ders_basis_funs<T>(vspan, dv, m_v_degree, v, m_v_knots_vector, nv);
				// Eigen::Matrix<T, point_size,  Eigen::Dynamic> temp;
				// temp.resize(point_size, m_v_degree + 1);
				const auto& control_points = m_surface->get_control_points_ref();
				for (int k = 0; k <= du; ++k)
				{
					temp.setConstant(0.0);
					for (int s = 0; s <= m_v_degree; ++s)
					{
						temp.col(s) = control_points[vspan - m_v_degree + s].block(0, uspan - m_u_degree, point_size, m_u_degree + 1) * nu.col(k);
					}
					int dd = std::min(n - k, dv);
					for (int l = 0; l <= dd; ++l)
					{
						ration_result(k, l) = temp * nv.col(l);
					}
				}
				result = std::move(project_derivs_point<T, is_rational, point_size, n>::project_point_to_euclidean_space(ration_result));
				return ENUM_NURBS::NURBS_SUCCESS;

			}
		}

        ENUM_NURBS point_on_surface(surface_type* surf, T u, T v, Eigen::Vector<T, dim> &point)
        {
            int uspan = -1, vspan = -1;
			int u_degree = surf->m_u_degree;
			int v_degree = surf->m_v_degree;
			const auto& u_knots = surf->m_u_knots_vector;
			const auto& v_knots = surf->m_v_knots_vector;
			find_span<T>(u, u_degree, u_knots, u_degree, u_knots.rows() - u_degree - 1, uspan);
			find_span<T>(v, v_degree, v_knots, v_degree, v_knots.rows() - v_degree - 1, vspan);
            basis_functions<T>(uspan, u, u_degree, u_knots, u_basis_funs_array);
            basis_functions<T>(vspan, v, v_degree, v_knots, v_basis_funs_array);

			const auto& control_points = surf->m_control_points;
            for (int index = 0; index <= v_degree; ++index)
            {
                temp.col(index) = control_points[vspan + index - v_degree].middleCols(uspan - u_degree, u_degree + 1) * u_basis_funs_array.head(u_degree + 1);
            }
            if constexpr (is_rational == false)
            {
				point = temp.block(0, 0, point_size, v_degree + 1)* v_basis_funs_array.head(v_degree + 1);
            }
            else
            {
				point = project_point<T, is_rational, point_size>::project_point_to_euclidean_space(temp.block(0, 0, point_size, v_degree + 1) * v_basis_funs_array);
            }

            return ENUM_NURBS::NURBS_SUCCESS;
        }
	};


	template<typename surface_type>
	struct surfs_patch
	{
		using T = typename surface_type::Type;
		static constexpr int dim = surface_type::dimension;
		static_assert(dim == 3, "dim != 3");

		surface_type* m_surf;
		Box<T, dim> m_u_tangent_box;
		Box<T, dim> m_v_tangent_box;
		Box<T, dim> m_normal_box;

		Box<T, 2> m_uv_box;
		Box<T, dim> m_box;
		std::vector<surfs_patch*> m_children;
		std::vector<surfs_patch*> m_int_patches;
		bool m_box_is_valid{ false };
		~surfs_patch()
		{
			if (m_surf != nullptr)
			{
				delete m_surf;
			}
		}
	};

	template<typename left_surf_type, typename right_surf_type>
	struct int_surfs_pair
	{
		using T = typename left_surf_type::Type;
		static constexpr int dim = left_surf_type::dimension;
		static_assert(dim == 3, "dim != 3");
		static_assert(std::is_same_v<T, typename right_surf_type::Type>, "std::is_same_v<T, typename right_surf_type::Type>");
		static_assert(dim == right_surf_type::dimension, "dim != right_surf_type::dimension");

		surfs_patch<left_surf_type>* m_left_surfs_patch;
		surfs_patch<right_surf_type>* m_right_surfs_patch;

		~int_surfs_pair()
		{
			std::vector<surfs_patch<left_surf_type>*> current_patch{ m_left_surfs_patch };
			std::vector<surfs_patch<left_surf_type>*> next_childrens;
			while (current_patch.size() > 0)
			{
				for (int index = 0; index < current_patch.size(); ++index)
				{
					next_childrens.insert(next_childrens.end(), current_patch[index]->m_children.begin(), current_patch[index]->m_children.end());
					delete current_patch[index];
				}
				std::swap(current_patch, next_childrens);
				next_childrens.clear();
			}
			std::vector<surfs_patch<right_surf_type>*> right_current_patch{ m_right_surfs_patch };
			std::vector<surfs_patch<right_surf_type>*> right_next_childrens;
			while (right_current_patch.size() > 0)
			{
				for (int index = 0; index < right_current_patch.size(); ++index)
				{
					right_next_childrens.insert(right_next_childrens.end(), right_current_patch[index]->m_children.begin(), right_current_patch[index]->m_children.end());
					delete right_current_patch[index];
				}
				std::swap(right_current_patch, right_next_childrens);
				right_next_childrens.clear();
			}
		}
	};



	//计算第一基本型和第二基本型
	template<typename T, int dim>
	ENUM_NURBS eval_first_and_second_fundenmental(const Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> ders, const Eigen::Vector<T, dim>& normal,
		Eigen::Matrix2<T>& first, Eigen::Matrix2<T>& second)
	{
		static_assert(3 == dim, "3 != dim");
		first(0, 0) = ders(1, 0).dot(ders(1, 0));
		first(1, 0) = first(0, 1) = ders(1, 0).dot(ders(0, 1));
		first(1, 1) = ders(0, 1).dot(ders(0, 1));

		second(0, 0) = ders(2, 0).dot(normal);
		second(1, 0) = second(0, 1) = ders(1, 1).dot(normal);
		second(1, 1) = ders(0, 2).dot(normal);
		return ENUM_NURBS::NURBS_SUCCESS;
	}
	
	template<typename T, int dim>
	ENUM_NURBS eval_first_fundenmental(const Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> ders, Eigen::Matrix2<T>& first)
	{
		static_assert(3 == dim, "3 != dim");
		first(0, 0) = ders(1, 0).dot(ders(1, 0));
		first(1, 0) = first(0, 1) = ders(1, 0).dot(ders(0, 1));
		first(1, 1) = ders(0, 1).dot(ders(0, 1));

		return ENUM_NURBS::NURBS_SUCCESS;
	}

	template<typename surface_type>
	struct bounding_box
	{
		using T = typename surface_type::Type;
		static constexpr int dim = surface_type::dimension;
		static ENUM_NURBS get_ders_bounding_box(int degree, const Eigen::MatrixX<surface_type*>& surf, const Box<T, 2>& interval, Eigen::MatrixX<Box<T, dim>>& box)
		{
			assert(degree < surf.cols() && degree < surf.rows());
			box.resize(degree + 1, degree + 1);
			if constexpr (surface_type::is_ratio == false)
			{
				for (int i = 0; i <= degree; ++i)
				{
					for (int j = 0; j <= degree - i; ++j)
					{
						surface_type sub_surface;
						Box<T, 2> interval_copy = interval;
						surf(i, j)->sub_divide(interval_copy, sub_surface);
						sub_surface.get_box(box(i, j));
					}
				}
				return ENUM_NURBS::NURBS_SUCCESS;
			}
			else
			{
				ENUM_NURBS error_code;
				Eigen::MatrixX<int> bin = binary_coeff(degree + 1);
				Eigen::MatrixX<Box<T, 1>> w_box;
				w_box.resize(degree + 1, degree + 1);
				for (int i = 0; i <= degree; ++i)
				{
					for (int j = 0; j <= degree - i; ++j)
					{
						surface_type sub_surface;
						Box<T, 2> interval_copy = interval;
						surf(i, j)->sub_divide(interval_copy, sub_surface);
						Box<T, dim + 1> A;
						sub_surface.get_ratio_box(A);
						Box<T, dim> Aij;
						Aij.Min = A.Min.template block<dim, 1>(0, 0);
						Aij.Max = A.Max.template block<dim, 1>(0, 0);
						w_box(i, j).Min[0] = A.Min[dim];
						w_box(i, j).Max[0] = A.Max[dim];
						for (int k = 1; k <= i; ++k)
						{
							Aij = Aij - box(i - k, j) * (w_box(k, 0) * bin(i, k));
						}
						for (int k = 1; k <= j; ++k)
						{
							Aij = Aij - box(i, j - k) * (w_box(0, k) * bin(j, k));
						}
						for (int k = 1; k <= i; ++k)
						{
							for (int r = 1; r <= j; ++r)
							{
								Aij = Aij - box(k, r) * (w_box(k, r) * (bin(i, k) * bin(j, r)));
							}
						}
						error_code = interval_algorithm::divide(Aij, w_box(0, 0), box(i, j));
						if (error_code != ENUM_NURBS::NURBS_SUCCESS)
						{
							return error_code;
						}
					}
				}
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		static ENUM_NURBS get_ration_bcurve_normal_and_ders_box(const Box<T, dim>& box, const surface_type& surf, const surface_type& der10, const surface_type& der01, Box<T, dim>& box10,
														Box<T, dim>& box01, Box<T, dim>& normal_box)
		{
			static_assert(surface_type::is_ratio == true, "surface_type::is_ratio == false");
			Box<T, 1> w_box;
			surf.get_w_box(w_box);
			Box<T, 1> wd_box;

			Box<T, dim + 1> ratio_box;
			der01.get_ratio_box(ratio_box);
			Box<T, dim> A;
			A.Min = ratio_box.Min.template block<dim, 1>(0, 0);
			A.Max = ratio_box.Max.template block<dim, 1>(0, 0);
			wd_box.Min[0] = ratio_box.Min[dim];
			wd_box.Max[0] = ratio_box.Max[dim];
			ENUM_NURBS error_code = interval_algorithm::divide(A - box * wd_box, w_box, box01);
			if (error_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				return error_code;
			}

			der10.get_ratio_box(ratio_box);
			A.Min = ratio_box.Min.template block<dim, 1>(0, 0);
			A.Max = ratio_box.Max.template block<dim, 1>(0, 0);
			wd_box.Min[0] = ratio_box.Min[dim];
			wd_box.Max[0] = ratio_box.Max[dim];
			error_code = interval_algorithm::divide(A - box * wd_box, w_box, box10);
			if (error_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				return error_code;
			}

			normal_box = interval_algorithm::cross(box10, box01);
			return ENUM_NURBS::NURBS_SUCCESS;
		}
	};

	enum POINT_TYPE
	{
		UNKNOWN = -1,
		BOUNDARY = 0,
		SINGULAR = 1,
		LOOP = 2
	};

	//intersect points(以后有空再改吧)
	template<typename T, int dim>
	struct surf_surf_intersect_point
	{
		Eigen::Vector<T, dim> m_point;
		Eigen::Vector4<T> m_uv;
		Box<T, 4> m_priori_enclosure;
		double m_previous_step = 100000;
	};

	template<typename T, int dim>
	struct surfs_int_points_chat
	{
		std::vector<surf_surf_intersect_point<T, dim>> m_inter_points;
		bool m_is_transversal{ true };
		bool m_is_positive_direction{ true };
		POINT_TYPE start_point_type{ POINT_TYPE::UNKNOWN };
		POINT_TYPE end_point_type{ POINT_TYPE::UNKNOWN };
	};


	//面面相交的结果
	template<typename T, int dim>
	struct surf_surf_int
	{
		std::vector<surfs_int_points_chat<T, dim>> m_int_chats;
		std::vector<surf_surf_intersect_point<T, dim>> m_isolate_points;

		//每个交点链表的交点的个数
		~surf_surf_int() { }
	};


	// TODO: 整理和完善代码
	template<typename left_surface_type, typename right_surface_type>
	class nurbs_surfaces_intersect
	{
		using T = typename left_surface_type::Type;
		static constexpr int dim = left_surface_type::dimension;

		using right_curve = typename right_surface_type::iso_curve_type;
		using left_curve = typename left_surface_type::iso_curve_type;
		static_assert(std::is_same<typename left_surface_type::Type, typename right_surface_type::Type>::value, "left surface type is not same with right surface type");
		static_assert(left_surface_type::dimension == right_surface_type::dimension, "left surface dimension not equal right surface dimensin");

	private:
		Eigen::Matrix3<left_surface_type*> m_left_ders;
		Eigen::Matrix3<right_surface_type*> m_right_ders;

		left_surface_type* m_left_normal_surface;
		right_surface_type* m_right_normal_surface;
		T m_angle_eps;
		Box<T, 4> m_product_box;

		int looop_count = 0;

		int surf_index = 0;

		std::vector<surfs_int_points_chat<T, dim>> m_boundary_chats;

		//奇异点必然为起始点或者终止点
		//奇异点的意思是在此点的切平面上存在至少两个方向v1, v2, 使得两个相交的平面在这两个方向上任意阶偏导数相等(弧长参数下), 奇异点必然相切
		std::vector<surfs_int_points_chat<T, dim>> m_singular_chats;
		std::vector<Box<T, 4>> m_singular_boxes;

		//闭曲线的追踪点
		std::vector<surfs_int_points_chat<T, dim>> m_loop_chats;

		////闭曲线的在参数域的切向量角度的积分，简单闭曲线的这个值应该为2pi或者-2pi
		//std::vector<T> m_loop_angles;
		Eigen::Matrix3<Eigen::Matrix2<left_surface_type>> m_current_left_ders;
		Eigen::Matrix3<Eigen::Matrix2<right_surface_type>> m_current_right_ders;
		Eigen::Matrix2<left_surface_type>  m_current_left_normal_surface;
		Eigen::Matrix2<right_surface_type> m_current_right_normal_surface;

		Eigen::Matrix2<Box<T, 2>> m_current_left_param_box;
		Eigen::Matrix2<Box<T, 2>> m_current_right_param_box;
		std::array<Box<T, 1>, 2> alpha_d;
		std::array<Box<T, 1>, 2> beta_d;

		Box<T, dim> intersect_box;

		surf_compute<left_surface_type>* m_left_surf_compute;
		surf_compute<right_surface_type>* m_right_surf_compute;


		using left_surface_type_sptr = std::shared_ptr<left_surface_type>;
		using right_surface_type_sptr = std::shared_ptr<right_surface_type>;


	public:
		surf_surf_int<T, dim> m_result;


#ifdef DSSI
		std::vector<Eigen::Vector3<T>> loop_points;
#endif // DSSI

		// debug used
		bool is_point_in_intcurve(const Eigen::Vector4<T>& param, bool is_transversal, size_t& curve_index)
		{
			Eigen::Vector2<T> test_param = param.template block<2, 1>(2, 0);
			curve_index = 0;
			for (surfs_int_points_chat<T, dim>& int_chat : m_result.m_int_chats)
			{
				if (is_transversal != int_chat.m_is_transversal)
				{
					continue;
				}
				int int_points_count = int_chat.m_inter_points.size();
				for (int index = 0; index < int_points_count - 1; ++index)
				{
					const Box<T, 4>& priori_enclosure = int_chat.m_inter_points[index].m_priori_enclosure;
					if (priori_enclosure.is_contain_point(param))
					{
						Eigen::Vector4<T> current_param = int_chat.m_inter_points[index].m_uv;
						Eigen::Vector2<T> next_param = int_chat.m_inter_points[index + 1].m_uv.template block<2, 1>(2, 0);
						int type = is_transversal == true ? 0 : 3;
						if (index == 0 && int_chat.start_point_type == POINT_TYPE::SINGULAR)
						{
							next_param = int_chat.m_inter_points[index].m_uv.template block<2, 1>(2, 0);
							current_param = int_chat.m_inter_points[index + 1].m_uv;

							if (is_point_on_arc(current_param, next_param, test_param, priori_enclosure, type, int_chat.m_is_positive_direction == false) == true)
							{
								return true;
							}
						}
						else
						{
							if (is_point_on_arc(current_param, next_param, test_param, priori_enclosure, type, int_chat.m_is_positive_direction) == true)
							{
								return true;
							}
						}
					}
				}
				curve_index += 1;
			}
			return false;
		}

	public:
		nurbs_surfaces_intersect() {};
		~nurbs_surfaces_intersect()
		{
			if (m_left_surf_compute != nullptr)
			{
				delete m_left_surf_compute;
			}
			if (m_right_surf_compute != nullptr)
			{
				delete m_right_surf_compute;
			}
			for (Eigen::Index col = 0; col < 3; ++col)
			{
				for (Eigen::Index row = 0; row < 3; ++row)
				{
					left_surface_type* left_surface = m_left_ders(row, col);
					if (left_surface != nullptr)
					{
						delete left_surface;
					}
					right_surface_type* right_surface = m_right_ders(row, col);
					if (right_surface != nullptr)
					{
						delete right_surface;
					}
				}
			}
		}

		ENUM_NURBS init(left_surface_type* left_surface, right_surface_type* right_surface, T angle_eps = 0.2)
		{
			m_left_ders(0, 0) = new left_surface_type(*left_surface);
			m_left_ders(1, 0) = new left_surface_type();
			m_left_ders(2, 0) = new left_surface_type();
			m_left_ders(0, 1) = new left_surface_type();
			m_left_ders(1, 1) = new left_surface_type();
			m_left_ders(2, 1) = nullptr;
			m_left_ders(0, 2) = new left_surface_type();
			m_left_ders(1, 2) = nullptr;
			m_left_ders(2, 2) = nullptr;

			m_right_ders(0, 0) = new right_surface_type(*right_surface);
			m_right_ders(1, 0) = new right_surface_type();
			m_right_ders(2, 0) = new right_surface_type();
			m_right_ders(0, 1) = new right_surface_type();
			m_right_ders(1, 1) = new right_surface_type();
			m_right_ders(2, 1) = nullptr;
			m_right_ders(0, 2) = new right_surface_type();
			m_right_ders(1, 2) = nullptr;
			m_right_ders(2, 2) = nullptr;

			m_left_ders(0, 0)->tangent_u_surface(*m_left_ders(1, 0));
			m_left_ders(0, 0)->tangent_v_surface(*m_left_ders(0, 1));

			//Xuu
			m_left_ders(1, 0)->tangent_u_surface(*(m_left_ders(2, 0)));
			//必须保证Xuv = Xvu, 对于beizer是必然满足的
			m_left_ders(1, 0)->tangent_v_surface(*(m_left_ders(1, 1)));
			//Xvv
			m_left_ders(0, 1)->tangent_v_surface(*(m_left_ders(0, 2)));

			m_right_ders(0, 0)->tangent_u_surface(*m_right_ders(1, 0));
			m_right_ders(0, 0)->tangent_v_surface(*m_right_ders(0, 1));

			//Yuu
			m_right_ders(1, 0)->tangent_u_surface(*(m_right_ders(2, 0)));
			//必须保证Yuv = Yvu, 对于beizer是必然满足的
			m_right_ders(1, 0)->tangent_v_surface(*(m_right_ders(1, 1)));
			//Yvv
			m_right_ders(0, 1)->tangent_v_surface(*(m_right_ders(0, 2)));

			if constexpr (left_surface_type::is_ratio == false)
			{
				m_left_normal_surface = new left_surface_type();
				m_left_ders(0, 0)->eval_normal_surface(*m_left_normal_surface);
			}
			if constexpr (right_surface_type::is_ratio == false)
			{
				m_right_normal_surface = new right_surface_type();
				m_right_ders(0, 0)->eval_normal_surface(*m_right_normal_surface);
			}

			m_angle_eps = angle_eps;
			m_product_box = left_surface->get_interval().product_box(right_surface->get_interval());

			m_left_surf_compute = new surf_compute<left_surface_type>(left_surface);
			m_right_surf_compute = new surf_compute<right_surface_type>(right_surface);


			return ENUM_NURBS::NURBS_SUCCESS;
		}

		// type = 0, 1, 2
		ENUM_NURBS eval_tangent(const Eigen::Vector<T, 4>& initial_param, Eigen::Vector<T, dim>& tangent,
			Eigen::Vector<T, 4>& param_tangent, int type)
		{
			Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> left_ders;
			m_left_surf_compute<1>(initial_param[0], initial_param[1], left_ders);
			// m_left_ders(0, 0)->template derivative_on_surface<1>(initial_param[0], initial_param[1], left_ders);
			Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> right_ders;
			m_right_surf_compute<1>(initial_param[2], initial_param[3], right_ders);
			// m_right_ders(0, 0)->template derivative_on_surface<1>(initial_param[2], initial_param[3], right_ders);

			Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
			T left_normal_len = left_normal.norm();
			left_normal.normalize();

			Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
			T right_normal_len = right_normal.norm();
			right_normal.normalize();
			if (type == 0)
			{
				tangent = left_normal.cross(right_normal);
				tangent.normalize();
			}
			else
			{
				Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_ders;
				Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_ders;
				{
					// 这里和上面有重复计算，之后优化
					m_left_surf_compute->template derivative_on_surface<2>(initial_param[0], initial_param[1], left_ders);
					m_right_surf_compute->template derivative_on_surface<2>(initial_param[2], initial_param[3], right_ders);
					
					// m_left_ders(0, 0)->template derivative_on_surface<2>(initial_param[0], initial_param[1], left_ders);
					// m_right_ders(0, 0)->template derivative_on_surface<2>(initial_param[2], initial_param[3], right_ders);
				}

				T denominator = 1.0 / right_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal);
				T a11 = left_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal) * denominator;
				T a12 = left_ders(0, 1).cross(right_ders(0, 1)).dot(right_normal) * denominator;
				T a21 = right_ders(1, 0).cross(left_ders(1, 0)).dot(right_normal) * denominator;
				T a22 = right_ders(1, 0).cross(left_ders(0, 1)).dot(right_normal) * denominator;
				T left_L = left_ders(2, 0).dot(left_normal);
				T left_M = left_ders(1, 1).dot(left_normal);
				T left_N = left_ders(0, 2).dot(left_normal);

				T right_L = right_ders(2, 0).dot(right_normal);
				T right_M = right_ders(1, 1).dot(right_normal);
				T right_N = right_ders(0, 2).dot(right_normal);

				T b11 = a11 * a11 * right_L;
				b11 += 2 * a11 * a21 * right_M;
				b11 += a21 * a21 * right_N;
				b11 -= left_L;

				T b12 = a11 * a12 * right_L;
				b12 += (a11 * a22 + a21 * a12) * right_M;
				b12 += a21 * a22 * right_N;
				b12 -= left_M;

				T b22 = a12 * a12 * right_L;
				b22 += 2 * a12 * a22 * right_M;
				b22 += a22 * a22 * right_N;
				b22 -= left_N;
				//double delta = b12 * b12 - b11 * b22;
				if (type == 1)
				{
					T nu = b12 / b11;
					tangent = left_ders(0, 1) - nu * left_ders(1, 0);
					tangent.normalize();
					//return ENUM_NURBS::NURBS_SUCCESS;
				}
				else
				{
					T mu = b12 / b22;
					tangent = left_ders(1, 0) - mu * left_ders(0, 1);
					tangent.normalize();
					//return ENUM_NURBS::NURBS_SUCCESS;
				}
			}

			param_tangent[0] = tangent.cross(left_ders(0, 1)).dot(left_normal) / left_normal_len;
			param_tangent[1] = left_ders(1, 0).cross(tangent).dot(left_normal) / left_normal_len;
			param_tangent[2] = tangent.cross(right_ders(0, 1)).dot(right_normal) / right_normal_len;
			param_tangent[3] = right_ders(1, 0).cross(tangent).dot(right_normal) / right_normal_len;
			return ENUM_NURBS::NURBS_SUCCESS;
		}
	
		
		template<typename surface_type>
		ENUM_NURBS eval_normal_and_ders_box(const surface_type& surface, Box<T, dim>& u_tangent_box, Box<T, dim>& v_tangent_box, Box<T, dim>& normal_box)
		{
			if constexpr (surface_type::is_ratio == false)
			{
				surface.tangent_v_surface_box(v_tangent_box);
				surface.tangent_u_surface_box(u_tangent_box);
				normal_box = interval_algorithm::cross(u_tangent_box, v_tangent_box);
				// m_left_normal_surface->get_box(left_param_box, left_normal_box);
			}
			else
			{
				assert(false);
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		}
		
		template<typename surface_type>
		ENUM_NURBS eval_normal_and_ders_box(const Box<T, dim>& box, surface_type& surf, surface_type& der10, surface_type& der01, Box<T, dim>& u_tangent_box, Box<T, dim>& v_tangent_box, Box<T, dim>& normal_box)
		{
			if constexpr (surface_type::is_ratio == false)
			{
				der10.get_box(u_tangent_box);
				der01.get_box(v_tangent_box);
				//使用精确的normal_surface和使用两个box叉乘得到的结果竟然基本一样？？？？？？？
				normal_box = interval_algorithm::cross(u_tangent_box, v_tangent_box);
				// m_left_normal_surface->get_box(left_param_box, left_normal_box);
			}
			else
			{
				ENUM_NURBS error_code;
				error_code = bounding_box<surface_type>::get_ration_bcurve_normal_and_ders_box(box, surf, der10, der01, u_tangent_box, v_tangent_box, normal_box);
				if (error_code != ENUM_NURBS::NURBS_SUCCESS)
				{
					return error_code;
				}
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		ENUM_NURBS eval_preiamge_and_space_ders(Box<T, dim>& left_u_tangent_box, Box<T, dim>& left_v_tangent_box,
			Box<T, dim>& right_u_tangent_box, Box<T, dim>& right_v_tangent_box, const Box<T, dim>& left_normal_box, Box<T, dim>& right_normal_box,
			Box<T, dim>& space_ders, Box<T, 4>& param_ders,  bool& may_arrived_singular)
		{
			Box<T, dim> c;
			c = interval_algorithm::cross(left_normal_box, right_normal_box);
			Box<T, 1> left_normal_len2 = interval_algorithm::dot(left_normal_box, left_normal_box);
			Box<T, 1> right_normal_len2 = interval_algorithm::dot(right_normal_box, right_normal_box);

			Box<T, dim> unit_c;
			ENUM_NURBS error_code = interval_algorithm::normalized(c, unit_c);
			if (error_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				may_arrived_singular = true;
				return error_code;
			}
			space_ders = c;

			// std::vector<Box<T, 1>> alpha_d(2);
			// std::vector<Box<T, 1>> beta_d(2);
			if (false == interval_algorithm::divide(interval_algorithm::dot(left_normal_box, interval_algorithm::cross(unit_c, left_v_tangent_box)), left_normal_len2, alpha_d[0]))
				return ENUM_NURBS::NURBS_ERROR;

			if (false == interval_algorithm::divide(interval_algorithm::dot(left_normal_box, interval_algorithm::cross(left_u_tangent_box, unit_c)), left_normal_len2, alpha_d[1]))
				return ENUM_NURBS::NURBS_ERROR;

			if (false == interval_algorithm::divide(interval_algorithm::dot(right_normal_box, interval_algorithm::cross(unit_c, right_v_tangent_box)), right_normal_len2, beta_d[0]))
				return ENUM_NURBS::NURBS_ERROR;

			if (false == interval_algorithm::divide(interval_algorithm::dot(right_normal_box, interval_algorithm::cross(right_u_tangent_box, unit_c)), right_normal_len2, beta_d[1]))
				return ENUM_NURBS::NURBS_ERROR;
			param_ders.set_index_interval(0, alpha_d[0]);
			param_ders.set_index_interval(1, alpha_d[1]);
			param_ders.set_index_interval(2, beta_d[0]);
			param_ders.set_index_interval(3, beta_d[1]);
			return ENUM_NURBS::NURBS_SUCCESS;
		}


		std::array<Box<T, dim>, 3>  get_box(const Box<T, 2>& domain, bool is_left_surface)
		{
			// intervals:
			// (0, 1), (1, 1)
			//      (sp)
			// (0, 0) (1, 0)
			std::array<Box<T, dim>, 3> result;
			bool result_is_valid{ false };
			Box<T, 2> int_domain;
			
			if (is_left_surface)
			{
				m_current_left_param_box(0, 0).intersect(domain, int_domain);
				
				left_surface_type sub_surf;
				//(0, 0)
				if (int_domain.Max[0] - int_domain.Min[0] > 100 * PRECISION<T>::value && int_domain.Max[1] - int_domain.Min[1] > 100 * PRECISION<T>::value)
				{
					result_is_valid = true;
					m_current_left_ders(1, 0)(0, 0).template split_at_param<3>(int_domain.Min, sub_surf);
					sub_surf.get_box(result[0]);
					
					m_current_left_ders(0, 1)(0, 0).template split_at_param<3>(int_domain.Min, sub_surf);
					sub_surf.get_box(result[1]);
					
					m_current_left_normal_surface(0, 0).template split_at_param<3>(int_domain.Min, sub_surf);
					sub_surf.get_box(result[2]);
				}

				// (0, 1)
				m_current_left_param_box(0, 1).intersect(domain, int_domain);
				if (int_domain.Max[0] - int_domain.Min[0] > 100 * PRECISION<T>::value && int_domain.Max[1] - int_domain.Min[1] > 100 * PRECISION<T>::value)
				{
					Box<T, dim> temp;
					Eigen::Vector2<T> param(int_domain.Max[0], int_domain.Min[1]);
					m_current_left_ders(1, 0)(0, 1).template split_at_param<2>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[0].unit(temp);
						// result[0] = result[0] + temp;
					}
					else
					{
						result[0] = temp;
					}
					
					m_current_left_ders(0, 1)(0, 1).template split_at_param<2>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[1].unit(temp);
					}
					else
					{
						result[1] = temp;
					}
					
					m_current_left_normal_surface(0, 1).template split_at_param<2>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[2].unit(temp);
					}
					else
					{
						result[2] = temp;
					}
					result_is_valid = true;
				}

				// (1, 0)
				m_current_left_param_box(1, 0).intersect(domain, int_domain);
				if (int_domain.Max[0] - int_domain.Min[0] > 100 * PRECISION<T>::value && int_domain.Max[1] - int_domain.Min[1] > 100 * PRECISION<T>::value)
				{
					Box<T, dim> temp;
					Eigen::Vector2<T> param(int_domain.Min[0], int_domain.Max[1]);
					m_current_left_ders(1, 0)(1, 0).template split_at_param<1>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[0].unit(temp);
					}
					else
					{
						result[0] = temp;
					}
					
					m_current_left_ders(0, 1)(1, 0).template split_at_param<1>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[1].unit(temp);
					}
					else
					{
						result[1] = temp;
					}
					
					m_current_left_normal_surface(1, 0).template split_at_param<1>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[2].unit(temp);
					}
					else
					{
						result[2] = temp;
					}
					result_is_valid = true;
				}
				
				// (1, 1)
				m_current_left_param_box(1, 1).intersect(domain, int_domain);
				if (int_domain.Max[0] - int_domain.Min[0] > 100 * PRECISION<T>::value && int_domain.Max[1] - int_domain.Min[1] > 100 * PRECISION<T>::value)
				{
					Box<T, dim> temp;
					Eigen::Vector2<T> param(int_domain.Max[0], int_domain.Max[1]);
					m_current_left_ders(1, 0)(1, 1).template split_at_param<0>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[0].unit(temp);
					}
					else
					{
						result[0] = temp;
					}
					
					m_current_left_ders(0, 1)(1, 1).template split_at_param<0>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[1].unit(temp);
					}
					else
					{
						result[1] = temp;
					}
					
					m_current_left_normal_surface(1, 1).template split_at_param<0>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[2].unit(temp);
					}
					else
					{
						result[2] = temp;
					}
					result_is_valid = true;
				}
			}
			else
			{
				m_current_right_param_box(0, 0).intersect(domain, int_domain);

				right_surface_type sub_surf;
				//(0, 0)
				if (int_domain.Max[0] - int_domain.Min[0] > 100 * PRECISION<T>::value && int_domain.Max[1] - int_domain.Min[1] > 100 * PRECISION<T>::value)
				{
					result_is_valid = true;
					m_current_right_ders(1, 0)(0, 0).template split_at_param<3>(int_domain.Min, sub_surf);
					sub_surf.get_box(result[0]);

					m_current_right_ders(0, 1)(0, 0).template split_at_param<3>(int_domain.Min, sub_surf);
					sub_surf.get_box(result[1]);

					m_current_right_normal_surface(0, 0).template split_at_param<3>(int_domain.Min, sub_surf);
					sub_surf.get_box(result[2]);
				}

				// (0, 1)
				m_current_right_param_box(0, 1).intersect(domain, int_domain);
				if (int_domain.Max[0] - int_domain.Min[0] > 100 * PRECISION<T>::value && int_domain.Max[1] - int_domain.Min[1] > 100 * PRECISION<T>::value)
				{
					Box<T, dim> temp;
					Eigen::Vector2<T> param(int_domain.Max[0], int_domain.Min[1]);
					m_current_right_ders(1, 0)(0, 1).template split_at_param<2>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[0].unit(temp);
					}
					else
					{
						result[0] = temp;
					}

					m_current_right_ders(0, 1)(0, 1).template split_at_param<2>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[1].unit(temp);
					}
					else
					{
						result[1] = temp;
					}

					m_current_right_normal_surface(0, 1).template split_at_param<2>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[2].unit(temp);
					}
					else
					{
						result[2] = temp;
					}
					result_is_valid = true;
				}

				// (1, 0)
				m_current_right_param_box(1, 0).intersect(domain, int_domain);
				if (int_domain.Max[0] - int_domain.Min[0] > 100 * PRECISION<T>::value && int_domain.Max[1] - int_domain.Min[1] > 100 * PRECISION<T>::value)
				{
					Box<T, dim> temp;
					Eigen::Vector2<T> param(int_domain.Min[0], int_domain.Max[1]);
					m_current_right_ders(1, 0)(1, 0).template split_at_param<1>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[0].unit(temp);
					}
					else
					{
						result[0] = temp;
					}

					m_current_right_ders(0, 1)(1, 0).template split_at_param<1>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[1].unit(temp);
					}
					else
					{
						result[1] = temp;
					}

					m_current_right_normal_surface(1, 0).template split_at_param<1>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[2].unit(temp);
					}
					else
					{
						result[2] = temp;
					}
					result_is_valid = true;
				}

				// (1, 1)
				m_current_right_param_box(1, 1).intersect(domain, int_domain);
				if (int_domain.Max[0] - int_domain.Min[0] > 100 * PRECISION<T>::value && int_domain.Max[1] - int_domain.Min[1] > 100 * PRECISION<T>::value)
				{
					Box<T, dim> temp;
					Eigen::Vector2<T> param(int_domain.Max[0], int_domain.Max[1]);
					m_current_right_ders(1, 0)(1, 1).template split_at_param<0>(param, sub_surf);

					// ENUM_NURBS sub_divide(Box<T, 2>& uv_box, nurbs_surface<T, dim, -1, -1, -1, -1, is_rational>& sub_nurbs) const
					// right_surface_type temp_surface;
					// m_right_ders(1, 0)->sub_divide(m_current_right_param_box(1, 1), temp_surface);


					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[0].unit(temp);
					}
					else
					{
						result[0] = temp;
					}

					m_current_right_ders(0, 1)(1, 1).template split_at_param<0>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[1].unit(temp);
					}
					else
					{
						result[1] = temp;
					}

					m_current_right_normal_surface(1, 1).template split_at_param<0>(param, sub_surf);
					sub_surf.get_box(temp);
					if (result_is_valid)
					{
						result[2].unit(temp);
					}
					else
					{
						result[2] = temp;
					}
					result_is_valid = true;
				}
			}
			assert(result_is_valid == true);
			return result;
		}


		//2阶微分的box太大，要么能够缩小2阶微分的box，要么只使用一阶微分的box
		//validated solution of initial value可以计算交线在此部分的包围盒,此包围盒应该会有其他的一些用处
		template<unsigned order>
		ENUM_NURBS eval_preiamge_and_space_ders(const Box<T, 4>& param_box, std::array<Box<T, 4>, order>& param_ders,
			std::array<Box<T, dim>, order>& space_ders, bool& may_arrived_singular)
		{
			//目前仅仅支持三维
			static_assert(3 == dim, "dim != 3 is not supported");
			// static_assert(order == 1 || order == 2, "order != 1 && order != 2");
			static_assert(order == 1, "order != 1");
			may_arrived_singular = false;
			//计算交线的切向, 两个曲面非相切的时候;
			Box<T, dim> left_normal_box, right_normal_box;
			Box<T, dim> left_u_tangent_box, left_v_tangent_box;
			Box<T, dim> right_u_tangent_box, right_v_tangent_box;

			Box<T, 2> left_param_box(param_box.Min.template block<2, 1>(0, 0), param_box.Max.template block<2, 1>(0, 0));
			Box<T, 2> right_param_box(param_box.Min.template block<2, 1>(2, 0), param_box.Max.template block<2, 1>(2, 0));
			if constexpr (left_surface_type::is_ratio == false)
			{
				looop_count += 1;
				m_left_surf_compute->get_2_ders_sub_box2(left_param_box, left_u_tangent_box, left_v_tangent_box);
				// m_left_ders(0, 0)->get_2_ders_sub_box(left_param_box, left_u_tangent_box, left_v_tangent_box);
				left_normal_box = interval_algorithm::cross(left_u_tangent_box, left_v_tangent_box);
				//使用精确的normal_surface和使用两个box叉乘得到的结果竟然基本一样？？？？？？？
				// m_left_normal_surface->get_box(left_param_box, left_normal_box);
			}
			else
			{
				Eigen::MatrixX<Box<T, dim>> left_boxes;
				ENUM_NURBS error_code;
				error_code = bounding_box<left_surface_type>::get_ders_bounding_box(1, m_left_ders, left_param_box, left_boxes);
				if (error_code != ENUM_NURBS::NURBS_SUCCESS)
				{
					return error_code;
				}
				left_u_tangent_box = left_boxes(1, 0);
				left_v_tangent_box = left_boxes(0, 1);
				left_normal_box = interval_algorithm::cross(left_u_tangent_box, left_v_tangent_box);
			}
			if constexpr (right_surface_type::is_ratio == false)
			{
				looop_count += 1;
				// m_right_ders(0, 0)->get_2_ders_sub_box(right_param_box, right_u_tangent_box, right_v_tangent_box);
				m_right_surf_compute->get_2_ders_sub_box2(right_param_box, right_u_tangent_box, right_v_tangent_box);

				right_normal_box = interval_algorithm::cross(right_u_tangent_box, right_v_tangent_box);
				//使用精确的normal_surface和使用两个box叉乘得到的结果竟然基本一样？？？？？？？
				// m_right_normal_surface->get_box(right_param_box, right_normal_box);
				// Eigen::Vector3<T> vec1 = boxes[0].Min - right_u_tangent_box.Min;
				// Eigen::Vector3<T> vec2 = boxes[0].Max - right_u_tangent_box.Max;
				// if (vec1.norm() > 1e-4 || vec2.norm() > 1e-4)
				// {
				// 	std::cout << "error";
				// }
				// vec1 = boxes[1].Min - right_v_tangent_box.Min;
				// vec2 = boxes[1].Max - right_v_tangent_box.Max;
				// if (vec1.norm() > 1e-4 || vec2.norm() > 1e-4)
				// {
				// 	std::cout << "error";
				// }
				// vec1 = boxes[2].Min - right_normal_box.Min;
				// vec2 = boxes[2].Max - right_normal_box.Max;
				// if (vec1.norm() > 1e-4 || vec2.norm() > 1e-4)
				// {
				// 	std::cout << "error";
				// }
				// int index = 0;
			}
			else
			{
				Eigen::MatrixX<Box<T, dim>> right_boxes;
				ENUM_NURBS error_code;
				error_code = bounding_box<right_surface_type>::get_ders_bounding_box(1, m_right_ders, right_param_box, right_boxes);
				if (error_code != ENUM_NURBS::NURBS_SUCCESS)
				{
					return error_code;
				}
				right_u_tangent_box = right_boxes(1, 0);
				right_v_tangent_box = right_boxes(0, 1);
				right_normal_box = interval_algorithm::cross(right_u_tangent_box, right_v_tangent_box);
			}


			Box<T, dim> c;
			c = interval_algorithm::cross(left_normal_box, right_normal_box);
			Box<T, 1> left_normal_len2 = interval_algorithm::dot(left_normal_box, left_normal_box);
			Box<T, 1> right_normal_len2 = interval_algorithm::dot(right_normal_box, right_normal_box);

			Box<T, dim> unit_c;
			ENUM_NURBS error_code = interval_algorithm::normalized(c, unit_c);
			if (error_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				may_arrived_singular = true;
				return error_code;
			}
			//space_ders.resize(1);
			space_ders[0] = c;
			//param_ders.resize(1);

			// std::vector<Box<T, 1>> alpha_d(2);
			// std::vector<Box<T, 1>> beta_d(2);
			if (false == interval_algorithm::divide(interval_algorithm::dot(left_normal_box, interval_algorithm::cross(unit_c, left_v_tangent_box)), left_normal_len2, alpha_d[0]))
				return ENUM_NURBS::NURBS_ERROR;

			if (false == interval_algorithm::divide(interval_algorithm::dot(left_normal_box, interval_algorithm::cross(left_u_tangent_box, unit_c)), left_normal_len2, alpha_d[1]))
				return ENUM_NURBS::NURBS_ERROR;

			if (false == interval_algorithm::divide(interval_algorithm::dot(right_normal_box, interval_algorithm::cross(unit_c, right_v_tangent_box)), right_normal_len2, beta_d[0]))
				return ENUM_NURBS::NURBS_ERROR;

			if (false == interval_algorithm::divide(interval_algorithm::dot(right_normal_box, interval_algorithm::cross(right_u_tangent_box, unit_c)), right_normal_len2, beta_d[1]))
				return ENUM_NURBS::NURBS_ERROR;
			param_ders[0].set_index_interval(0, alpha_d[0]);
			param_ders[0].set_index_interval(1, alpha_d[1]);
			param_ders[0].set_index_interval(2, beta_d[0]);
			param_ders[0].set_index_interval(3, beta_d[1]);
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		//相切的情形的一阶微分的box
		//validated solution of initial value可以计算交线在此部分的包围盒,此包围盒应该会有其他的一些用处
		//目前只能使用两个surface的box内部normal的变化作为一个步长选择的标准
		template<unsigned order>
		ENUM_NURBS eval_preiamge_and_space_ders_tangent(const Box<T, 4>& param_box, std::array<Box<T, 4>, order>& param_ders,
			std::array<Box<T, dim>, order>& space_ders, int& type)
		{
			//目前仅仅支持三维
			static_assert(3 == dim, "dim != 3 is not supported");
			static_assert(order == 1, "order != 1");
			Box<T, dim> left_normal_box, right_normal_box;
			Box<T, dim> left_u_tangent_box, left_v_tangent_box;
			Box<T, dim> right_u_tangent_box, right_v_tangent_box;
			Box<T, dim> Xuu_box, Xuv_box, Xvv_box, Yuu_box, Yuv_box, Yvv_box;
			Box<T, 2> left_param_box(param_box.Min.template block<2, 1>(0, 0), param_box.Max.template block<2, 1>(0, 0));
			Box<T, 2> right_param_box(param_box.Min.template block<2, 1>(2, 0), param_box.Max.template block<2, 1>(2, 0));
			if constexpr (left_surface_type::is_ratio == false)
			{
				m_left_ders(1, 0)->get_box(left_param_box, left_u_tangent_box);
				m_left_ders(0, 1)->get_box(left_param_box, left_v_tangent_box);
				m_left_ders(2, 0)->get_box(left_param_box, Xuu_box);
				m_left_ders(1, 1)->get_box(left_param_box, Xuv_box);
				m_left_ders(0, 2)->get_box(left_param_box, Xvv_box);

				//使用精确的normal_surface和使用两个box叉乘得到的结果竟然基本一样？？？？？？？
				m_left_normal_surface->get_box(left_param_box, left_normal_box);
			}
			else
			{
				Eigen::MatrixX<Box<T, dim>> left_boxes;
				ENUM_NURBS error_code;
				error_code = bounding_box<left_surface_type>::get_ders_bounding_box(2, m_left_ders, left_param_box, left_boxes);
				if (error_code != ENUM_NURBS::NURBS_SUCCESS)
				{
					return error_code;
				}
				left_u_tangent_box = left_boxes(1, 0);
				left_v_tangent_box = left_boxes(0, 1);
				Xuu_box = left_boxes(2, 0);
				Xuv_box = left_boxes(1, 1);
				Xvv_box = left_boxes(0, 2);

				left_normal_box = interval_algorithm::cross(left_u_tangent_box, left_v_tangent_box);
			}
			if constexpr (right_surface_type::is_ratio == false)
			{
				m_right_ders(1, 0)->get_box(right_param_box, right_u_tangent_box);
				m_right_ders(0, 1)->get_box(right_param_box, right_v_tangent_box);
				m_right_ders(2, 0)->get_box(right_param_box, Yuu_box);
				m_right_ders(1, 1)->get_box(right_param_box, Yuv_box);
				m_right_ders(0, 2)->get_box(right_param_box, Yvv_box);

				//使用精确的normal_surface和使用两个box叉乘得到的结果竟然基本一样？？？？？？？
				m_right_normal_surface->get_box(right_param_box, right_normal_box);
			}
			else
			{
				Eigen::MatrixX<Box<T, dim>> right_boxes;
				ENUM_NURBS error_code;
				error_code = bounding_box<right_surface_type>::get_ders_bounding_box(2, m_right_ders, right_param_box, right_boxes);
				if (error_code != ENUM_NURBS::NURBS_SUCCESS)
				{
					return error_code;
				}
				right_u_tangent_box = right_boxes(1, 0);
				right_v_tangent_box = right_boxes(0, 1);

				Yuu_box = right_boxes(2, 0);
				Yuv_box = right_boxes(1, 1);
				Yvv_box = right_boxes(0, 2);

				right_normal_box = interval_algorithm::cross(right_u_tangent_box, right_v_tangent_box);
			}

			// eval_point_interval(m_left_normal_surface, left_param_box, left_normal_box);
			// eval_point_interval(m_right_normal_surface, right_param_box, right_normal_box);

			Box<T, 1> left_normal_len2 = interval_algorithm::dot(left_normal_box, left_normal_box);
			Box<T, 1> right_normal_len2 = interval_algorithm::dot(right_normal_box, right_normal_box);

			Box<T, 1> normal_dot = interval_algorithm::dot(left_normal_box, right_normal_box);
			Eigen::Vector<T, 1> zero_num;
			zero_num[0] = 0;
			if (normal_dot.is_contain_point(zero_num, TDEFAULT_ERROR<T>::value) == true)
			{
				return ENUM_NURBS::NURBS_ERROR;
			}
			//int flag = normal_dot.Min[0] > 0 ? 1 : -1;
			Box<T, dim> conormal;
			ENUM_NURBS code = interval_algorithm::normalized(right_normal_box, conormal);
			if (ENUM_NURBS::NURBS_SUCCESS != code)
			{
				return code;
			}

			Eigen::Vector<T, dim> origin;
			origin.setConstant(0.0);
			T right_normal_min_length = right_normal_box.eval_minimal_distance(origin);
			T right_normal_max_length = right_normal_box.eval_minimal_distance(origin);
			Box<T, 1> denominator;
			denominator.Min[0] = right_normal_min_length;
			denominator.Max[0] = right_normal_max_length;
			Box<T, 1> denominator_r, unitInterval;
			unitInterval.Min[0] = 1.0;
			unitInterval.Max[0] = 1.0;
			code = interval_algorithm::divide(unitInterval, denominator, denominator_r);
			if (ENUM_NURBS::NURBS_SUCCESS != code)
			{
				return code;
			}

			Box<T, 1> a11 = interval_algorithm::dot(interval_algorithm::cross(left_u_tangent_box, right_v_tangent_box), conormal) * denominator_r;
			Box<T, 1> a12 = interval_algorithm::dot(interval_algorithm::cross(left_v_tangent_box, right_v_tangent_box), conormal) * denominator_r;
			Box<T, 1> a21 = interval_algorithm::dot(interval_algorithm::cross(right_u_tangent_box, left_u_tangent_box), conormal) * denominator_r;
			Box<T, 1> a22 = interval_algorithm::dot(interval_algorithm::cross(right_u_tangent_box, left_v_tangent_box), conormal) * denominator_r;

			Box<T, 1> left_L = interval_algorithm::dot(Xuu_box, conormal);
			Box<T, 1> left_M = interval_algorithm::dot(Xuv_box, conormal);
			Box<T, 1> left_N = interval_algorithm::dot(Xvv_box, conormal);

			Box<T, 1> right_L = interval_algorithm::dot(Yuu_box, conormal);
			Box<T, 1> right_M = interval_algorithm::dot(Yuv_box, conormal);
			Box<T, 1> right_N = interval_algorithm::dot(Yvv_box, conormal);

			Box<T, 1> b11 = a11 * a11;
			b11.Min[0] = std::max(0.0, b11.Min[0]);
			b11 = b11 * right_L;
			b11 = b11 + (a11 * a21 * right_M) * 2.0;
			Box<T, 1> temp = a21 * a21;
			temp.Min[0] = std::max(0.0, temp.Min[0]);
			b11 = b11 + temp * right_N;
			b11 = b11 - left_L;

			Box<T, 1> b12 = a11 * a12 * right_L;
			b12 = b12 + (a11 * a22 + a21 * a12) * right_M;
			b12 = b12 + a21 * a22 * right_N;
			b12 = b12 - left_M;

			Box<T, 1> b22 = a12 * a12;
			b22.Min[0] = std::max(0.0, b22.Min[0]);
			b22 = b22 * right_L;
			b22 = b22 + (a12 * a22 * right_M) * 2.0;
			temp = a22 * a22;
			temp.Min[0] = std::max(0.0, temp.Min[0]);
			b22 = b22 + temp * right_N;
			b22 = b22 - left_N;

			bool flag1 = b11.is_contain_point(zero_num, TDEFAULT_ERROR<T>::value);
			//bool flag2 = b12.is_contain_point(zero_num, TDEFAULT_ERROR<T>::value);
			bool flag2 = b22.is_contain_point(zero_num, TDEFAULT_ERROR<T>::value);
			if (flag1 && flag2)
			{
				return ENUM_NURBS::NURBS_ERROR;
			}
			T b11_extream = std::min(std::abs(b11.Min[0]), std::abs(b11.Max[0]));
			T b22_extream = std::min(std::abs(b22.Min[0]), std::abs(b22.Max[0]));

			if (flag2 == true || (b11_extream > b22_extream && flag1 == false))
			{
				Box<T, 1> nu;
				interval_algorithm::divide(b12, b11, nu);
				Box<T, dim> tangent = left_u_tangent_box * nu + left_v_tangent_box;
				code = interval_algorithm::normalized(tangent, space_ders[0]);
				if (ENUM_NURBS::NURBS_SUCCESS != code)
				{
					return code;
				}
				type = 1;
			}
			else
			{
				Box<T, 1> mu;
				interval_algorithm::divide(b12, b22, mu);
				Box<T, dim> tangent = left_v_tangent_box * mu + left_u_tangent_box;
				code = interval_algorithm::normalized(tangent, space_ders[0]);
				if (ENUM_NURBS::NURBS_SUCCESS != code)
				{
					return code;
				}
				type = 2;
			}
			std::vector<Box<T, 1>> alpha_d(2);
			std::vector<Box<T, 1>> beta_d(2);
			if (false == interval_algorithm::divide(interval_algorithm::dot(left_normal_box, interval_algorithm::cross(space_ders[0], left_v_tangent_box)), left_normal_len2, alpha_d[0]))
				return ENUM_NURBS::NURBS_ERROR;

			if (false == interval_algorithm::divide(interval_algorithm::dot(left_normal_box, interval_algorithm::cross(left_u_tangent_box, space_ders[0])), left_normal_len2, alpha_d[1]))
				return ENUM_NURBS::NURBS_ERROR;

			if (false == interval_algorithm::divide(interval_algorithm::dot(right_normal_box, interval_algorithm::cross(space_ders[0], right_v_tangent_box)), right_normal_len2, beta_d[0]))
				return ENUM_NURBS::NURBS_ERROR;

			if (false == interval_algorithm::divide(interval_algorithm::dot(right_normal_box, interval_algorithm::cross(right_u_tangent_box, space_ders[0])), right_normal_len2, beta_d[1]))
				return ENUM_NURBS::NURBS_ERROR;
			param_ders[0].set_index_interval(0, alpha_d[0]);
			param_ders[0].set_index_interval(1, alpha_d[1]);
			param_ders[0].set_index_interval(2, beta_d[0]);
			param_ders[0].set_index_interval(3, beta_d[1]);

			return ENUM_NURBS::NURBS_SUCCESS;
		}

		ENUM_NURBS eval_preiamge_and_space_ders2(const Eigen::Vector<T, 4>& param, std::vector<Eigen::Vector<T, 4>>& param_ders, std::vector<Eigen::Vector<T, dim>>& space_ders, int type)
		{
			static_assert(3 == dim, "3 != dim");

			param_ders.clear();
			space_ders.clear();
			Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_ders;
			Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_ders;
			if (type == 0)
			{
				// Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> left_ders_temp;
				// // m_left_ders(0, 0)->template derivative_on_surface<1, 1>(param[0], param[1], left_ders_temp);
				// m_left_surf_compute->template derivative_on_surface<1>(param[0], param[1], left_ders_temp);
				// for (int i = 0; i < 2; ++i)
				// {
				// 	for (int j = 0; j < 2; ++j)
				// 	{
				// 		left_ders(i, j) = left_ders_temp(i, j);
				// 	}
				// }

				// Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> right_ders_temp;
				// m_right_surf_compute->template derivative_on_surface<1>(param[2], param[3], right_ders_temp);
				// // m_right_ders(0, 0)->template derivative_on_surface<1, 1>(param[2], param[3], right_ders_temp);
				// for (int i = 0; i < 2; ++i)
				// {
				// 	for (int j = 0; j < 2; ++j)
				// 	{
				// 		right_ders(i, j) = right_ders_temp(i, j);
				// 	}
				// }
				if constexpr (left_surface_type::is_ratio == false && right_surface_type::is_ratio == false)
				{
					// m_left_ders(1, 0)->point_on_surface(param[0], param[1], left_ders(1, 0));
					m_left_surf_compute->point_on_surface(m_left_ders(1, 0), param[0], param[1], left_ders(1, 0));
					m_left_surf_compute->point_on_surface(m_left_ders(0, 1), param[0], param[1], left_ders(0, 1));
					m_right_surf_compute->point_on_surface(m_right_ders(1, 0), param[2], param[3], right_ders(1, 0));
					m_right_surf_compute->point_on_surface(m_right_ders(0, 1), param[2], param[3], right_ders(0, 1));
					// m_left_ders(1, 0)->point_on_surface(param[0], param[1], left_ders(1, 0));
					// m_left_ders(0, 1)->point_on_surface(param[0], param[1], left_ders(0, 1));
					// m_right_ders(1, 0)->point_on_surface(param[2], param[3], right_ders(1, 0));
					// m_right_ders(0, 1)->point_on_surface(param[2], param[3], right_ders(0, 1));
				}
				else
				{
					Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> left_ders_temp;
					// m_left_ders(0, 0)->template derivative_on_surface<1, 1>(param[0], param[1], left_ders_temp);
					m_left_surf_compute->template derivative_on_surface<1>(param[0], param[1], left_ders_temp);
					for (int i = 0; i < 2; ++i)
					{
						for (int j = 0; j < 2; ++j)
						{
							left_ders(i, j) = left_ders_temp(i, j);
						}
					}

					Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> right_ders_temp;
					m_right_surf_compute->template derivative_on_surface<1>(param[2], param[3], right_ders_temp);
					// m_right_ders(0, 0)->template derivative_on_surface<1, 1>(param[2], param[3], right_ders_temp);
					for (int i = 0; i < 2; ++i)
					{
						for (int j = 0; j < 2; ++j)
						{
							right_ders(i, j) = right_ders_temp(i, j);
						}
					}
				}
			}
			else
			{
				// m_left_ders(0, 0)->template derivative_on_surface<2, 2>(param[0], param[1], left_ders);
				// m_right_ders(0, 0)->template derivative_on_surface<2, 2>(param[2], param[3], right_ders);
				m_left_surf_compute->template derivative_on_surface<2>(param[0], param[1], left_ders);
				m_right_surf_compute->template derivative_on_surface<2>(param[2], param[3], right_ders);
			}

			Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
			left_normal.normalize();
			Eigen::Matrix2<T> left_first_fundenmental;
			eval_first_fundenmental(left_ders, left_first_fundenmental);

			Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
			right_normal.normalize();
			Eigen::Matrix2<T> right_first_fundenmental;
			eval_first_fundenmental(right_ders, right_first_fundenmental);

			Eigen::Vector<T, dim> tangent;
			if (type == 0)
			{
				tangent = left_normal.cross(right_normal);
				Eigen::Matrix<T, dim, 1> ders;
				ders.col(0) = tangent.normalized();
				space_ders.push_back(ders);
			}
			else
			{
				//int flag = right_normal.dot(left_normal) > 0 ? 1 : -1;
				T denominator = 1.0 / right_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal);
				T a11 = left_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal) * denominator;
				T a12 = left_ders(0, 1).cross(right_ders(0, 1)).dot(right_normal) * denominator;
				T a21 = right_ders(1, 0).cross(left_ders(1, 0)).dot(right_normal) * denominator;
				T a22 = right_ders(1, 0).cross(left_ders(0, 1)).dot(right_normal) * denominator;
				T left_L = left_ders(2, 0).dot(right_normal);
				T left_M = left_ders(1, 1).dot(right_normal);
				T left_N = left_ders(0, 2).dot(right_normal);

				T right_L = right_ders(2, 0).dot(right_normal);
				T right_M = right_ders(1, 1).dot(right_normal);
				T right_N = right_ders(0, 2).dot(right_normal);

				T b11 = a11 * a11 * right_L;
				b11 += 2 * a11 * a21 * right_M;
				b11 += a21 * a21 * right_N;
				b11 -= left_L;

				T b12 = a11 * a12 * right_L;
				b12 += (a11 * a22 + a21 * a12) * right_M;
				b12 += a21 * a22 * right_N;
				b12 -= left_M;

				T b22 = a12 * a12 * right_L;
				b22 += 2 * a12 * a22 * right_M;
				b22 += a22 * a22 * right_N;
				b22 -= left_N;
				double delta = 4.0 * (b12 * b12 - b11 * b22);
				bool is_branch = false;
				if (type == 3)
				{
					// 下面delta的值需要修改
					if (delta < -TDEFAULT_ERROR<T>::value)
					{
						T try_value = -b12 / b11;
						T value = b11 * (try_value * try_value) + 2 * b12 * try_value + b22;
						if (std::abs(value) < TDEFAULT_ERROR<T>::value)
						{
							delta = 0.0;
						}
						else
						{
							return ENUM_NURBS::NURBS_ISOLATED_TANGENTIAL_POINT;
						}
					}
					if (delta > TDEFAULT_ERROR<T>::value)
					{
						is_branch = true;
					}
					if (std::abs(b11) < TDEFAULT_ERROR<T>::value && std::abs(b12) < TDEFAULT_ERROR<T>::value
						&& std::abs(b22) < TDEFAULT_ERROR<T>::value)
					{
						return ENUM_NURBS::NURBS_HIGH_ORDER_TANGENTIAL;
					}
				}

				if (type == 1 || (type == 3 && std::abs(b11) >= std::abs(b22)))
				{
					std::vector<T> coeff;
					if (is_branch == true)
					{
						T delta_r = std::sqrt(delta);
						coeff.push_back(2 * b12 + delta_r);
						coeff.push_back(2 * b12 - delta_r);
					}
					else
					{
						coeff.push_back(2 * b12);
					}
					for (const T& coef : coeff)
					{
						T nu = coef / (2 * b11);
						tangent = left_ders(0, 1) - nu * left_ders(1, 0);
						Eigen::Matrix<T, dim, 1> ders;
						ders.col(0) = tangent.normalized();
						space_ders.push_back(ders);
					}
				}
				if (type == 2 || (type == 3 && std::abs(b11) < std::abs(b22)))
				{
					std::vector<T> coeff;
					if (is_branch == true)
					{
						T delta_r = std::sqrt(delta);
						coeff.push_back(2 * b12 + delta_r);
						coeff.push_back(2 * b12 - delta_r);
					}
					else
					{
						coeff.push_back(2 * b12);
					}
					for (const T& coef : coeff)
					{
						T mu = coef / (2 * b22);
						tangent = left_ders(1, 0) - mu * left_ders(0, 1);
						Eigen::Matrix<T, dim, 1> ders;
						ders.col(0) = tangent.normalized();
						space_ders.push_back(ders);
					}
					//return ENUM_NURBS::NURBS_SUCCESS;
				}
			}

			Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(left_first_fundenmental);
			Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd2(right_first_fundenmental);
			Eigen::Matrix2<T> mat;
			mat(0, 0) = mat(1, 1) = 1.0;
			mat(0, 1) = mat(1, 0) = left_normal.dot(right_normal);
			Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd3(mat);

			for (Eigen::Matrix<T, dim, 1>& space_der : space_ders)
			{
				Eigen::Vector<T, dim> tangent_vec = space_der.col(0);
				Eigen::Vector<T, 2> dot_ders;
				dot_ders[0] = tangent_vec.dot(left_ders(1, 0));
				dot_ders[1] = tangent_vec.dot(left_ders(0, 1));

				Eigen::Matrix<T, 4, 1> param_der;
				param_der.template block<2, 1>(0, 0) = matSvd.solve(dot_ders);
				if (matSvd.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;
				dot_ders[0] = tangent_vec.dot(right_ders(1, 0));
				dot_ders[1] = tangent_vec.dot(right_ders(0, 1));

				param_der.template block<2, 1>(2, 0) = matSvd2.solve(dot_ders);
				if (matSvd2.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;
				//计算两曲面的交线的二阶微分(弧长参数下)

				param_ders.push_back(param_der);
			}

			return ENUM_NURBS::NURBS_SUCCESS;
		}

		/// @brief 计算两个曲面交线在各自参数域的原像曲线的1-2阶微分
		/// @tparam surface_type 曲面类型
		/// @tparam T double, flaot...
		/// @tparam dim 维数(目前仅仅支持三维)
		/// @param left_surface 第一个曲面
		/// @param right_surface 第二个曲面
		/// @param u 第一个曲面的参数u
		/// @param v 第一个曲面的参数v
		/// @param s 第一个曲面的参数s
		/// @param t 第一个曲面的参数t
		/// @param param_ders param_ders.col(0)为原像曲线的一阶微分, param_ders.col(1)为原像曲线的二阶微分
		/// @param space_ders space_ders.col(0)为交线的点, space_ders.col(1)为交线的一阶微分, space_ders.col(1)为交先的二阶微分
		/// @return 错误码
		template<int order>
		ENUM_NURBS eval_preiamge_and_space_ders(const Eigen::Vector<T, 4>& param, std::vector<Eigen::Matrix<T, 4, order>>& param_ders, std::vector<Eigen::Matrix<T, dim, order>>& space_ders, 
												std::array<T, 4>& ks, std::array<T, 2>& left_EG, std::array<T, 2>& right_EG, int type)
		{
			static_assert(3 == dim, "3 != dim");

			param_ders.clear();
			space_ders.clear();
			Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_ders;
			// m_left_ders(0, 0)->template derivative_on_surface<2, 2>(param[0], param[1], left_ders);
			m_left_surf_compute->template derivative_on_surface<2>(param[0], param[1], left_ders);
			Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
			left_normal.normalize();
			Eigen::Matrix2<T> left_first_fundenmental;
			Eigen::Matrix2<T> left_second_fundenmental;
			eval_first_and_second_fundenmental(left_ders, left_normal, left_first_fundenmental, left_second_fundenmental);
			ks[0] = left_second_fundenmental(0, 0) / left_first_fundenmental(0, 0);
			ks[1] = left_second_fundenmental(1, 1) / left_first_fundenmental(1, 1);
			left_EG[0] = left_first_fundenmental(0, 0);
			left_EG[1] = left_first_fundenmental(1, 1);
			Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_ders;
			// m_right_ders(0, 0)->template derivative_on_surface<2, 2>(param[2], param[3], right_ders);
			m_right_surf_compute->template derivative_on_surface<2>(param[2], param[3], right_ders);
			Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
			right_normal.normalize();
			Eigen::Matrix2<T> right_first_fundenmental;
			Eigen::Matrix2<T> right_second_fundenmental;
			eval_first_and_second_fundenmental(right_ders, right_normal, right_first_fundenmental, right_second_fundenmental);
			
			ks[2] = right_second_fundenmental(0, 0) / right_first_fundenmental(0, 0);
			ks[3] = right_second_fundenmental(1, 1) / right_first_fundenmental(1, 1);
			right_EG[0] = right_first_fundenmental(0, 0);
			right_EG[1] = right_first_fundenmental(1, 1);
			Eigen::Vector<T, dim> tangent;
			if (type == 0)
			{
				tangent = left_normal.cross(right_normal);
				Eigen::Matrix<T, dim, order> ders;
				ders.col(0) = tangent.normalized();
				space_ders.push_back(ders);
			}
			else
			{
				//int flag = right_normal.dot(left_normal) > 0 ? 1 : -1;
				T denominator = 1.0 / right_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal);
				T a11 = left_ders(1, 0).cross(right_ders(0, 1)).dot(right_normal) * denominator;
				T a12 = left_ders(0, 1).cross(right_ders(0, 1)).dot(right_normal) * denominator;
				T a21 = right_ders(1, 0).cross(left_ders(1, 0)).dot(right_normal) * denominator;
				T a22 = right_ders(1, 0).cross(left_ders(0, 1)).dot(right_normal) * denominator;
				T left_L = left_ders(2, 0).dot(right_normal);
				T left_M = left_ders(1, 1).dot(right_normal);
				T left_N = left_ders(0, 2).dot(right_normal);

				T right_L = right_ders(2, 0).dot(right_normal);
				T right_M = right_ders(1, 1).dot(right_normal);
				T right_N = right_ders(0, 2).dot(right_normal);

				T b11 = a11 * a11 * right_L;
				b11 += 2 * a11 * a21 * right_M;
				b11 += a21 * a21 * right_N;
				b11 -= left_L;

				T b12 = a11 * a12 * right_L;
				b12 += (a11 * a22 + a21 * a12) * right_M;
				b12 += a21 * a22 * right_N;
				b12 -= left_M;

				T b22 = a12 * a12 * right_L;
				b22 += 2 * a12 * a22 * right_M;
				b22 += a22 * a22 * right_N;
				b22 -= left_N;
				double delta = 4.0 * (b12 * b12 - b11 * b22);
				bool is_branch = false;
				if (type == 3)
				{
					// 下面delta的值需要修改
					if (delta < -TDEFAULT_ERROR<T>::value)
					{
						T try_value = -b12 / b11;
						T value = b11 * (try_value * try_value) + 2 * b12 * try_value + b22;
						if (std::abs(value) < TDEFAULT_ERROR<T>::value)
						{
							delta = 0.0;
						}
						else
						{
							return ENUM_NURBS::NURBS_ISOLATED_TANGENTIAL_POINT;
						}
					}
					if (delta > TDEFAULT_ERROR<T>::value)
					{
						is_branch = true;
					}
					if (std::abs(b11) < TDEFAULT_ERROR<T>::value && std::abs(b12) < TDEFAULT_ERROR<T>::value
						&& std::abs(b22) < TDEFAULT_ERROR<T>::value)
					{
						return ENUM_NURBS::NURBS_HIGH_ORDER_TANGENTIAL;
					}
				}

				if (type == 1 || (type == 3 && std::abs(b11) >= std::abs(b22)))
				{
					std::vector<T> coeff;
					if (is_branch == true)
					{
						T delta_r = std::sqrt(delta);
						coeff.push_back(2 * b12 + delta_r);
						coeff.push_back(2 * b12 - delta_r);
					}
					else
					{
						coeff.push_back(2 * b12);
					}
					for (const T& coef : coeff)
					{
						T nu = coef / (2 * b11);
						tangent = left_ders(0, 1) - nu * left_ders(1, 0);
						Eigen::Matrix<T, dim, order> ders;
						ders.col(0) = tangent.normalized();
						space_ders.push_back(ders);
					}
				}
				if (type == 2 || (type == 3 && std::abs(b11) < std::abs(b22)))
				{
					std::vector<T> coeff;
					if (is_branch == true)
					{
						T delta_r = std::sqrt(delta);
						coeff.push_back(2 * b12 + delta_r);
						coeff.push_back(2 * b12 - delta_r);
					}
					else
					{
						coeff.push_back(2 * b12);
					}
					for (const T& coef : coeff)
					{
						T mu = coef / (2 * b22);
						tangent = left_ders(1, 0) - mu * left_ders(0, 1);
						Eigen::Matrix<T, dim, order> ders;
						ders.col(0) = tangent.normalized();
						space_ders.push_back(ders);
					}
					//return ENUM_NURBS::NURBS_SUCCESS;
				}
			}

			Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(left_first_fundenmental);
			Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd2(right_first_fundenmental);
			Eigen::Matrix2<T> mat;
			mat(0, 0) = mat(1, 1) = 1.0;
			mat(0, 1) = mat(1, 0) = left_normal.dot(right_normal);
			Eigen::JacobiSVD<Eigen::Matrix2<T>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd3(mat);

			for (Eigen::Matrix<T, dim, order>& space_der : space_ders)
			{
				Eigen::Vector<T, dim> tangent_vec = space_der.col(0);
				Eigen::Vector<T, 2> dot_ders;
				dot_ders[0] = tangent_vec.dot(left_ders(1, 0));
				dot_ders[1] = tangent_vec.dot(left_ders(0, 1));

				Eigen::Matrix<T, 4, order> param_der;
				param_der.template block<2, 1>(0, 0) = matSvd.solve(dot_ders);
				if (matSvd.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;
				dot_ders[0] = tangent_vec.dot(right_ders(1, 0));
				dot_ders[1] = tangent_vec.dot(right_ders(0, 1));

				param_der.template block<2, 1>(2, 0) = matSvd2.solve(dot_ders);
				if (matSvd2.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;
				//计算两曲面的交线的二阶微分(弧长参数下)

				if constexpr (order == 2)
				{
					//1.计算 d(u, v) 的第二基本型的值
					Eigen::Vector<T, 2> second_fundenmental_value;
					second_fundenmental_value[0] = param_der.template block<2, 1>(0, 0).col(0).transpose() * left_second_fundenmental * param_der.template block<2, 1>(0, 0).col(0);
					second_fundenmental_value[1] = param_der.template block<2, 1>(2, 0).transpose() * right_second_fundenmental * param_der.template block<2, 1>(2, 0);

					Eigen::Vector2<T> coeff = matSvd3.solve(second_fundenmental_value);
					if (matSvd3.info() != Eigen::Success)
						return ENUM_NURBS::NURBS_ERROR;

					//交线的二阶导数
					Eigen::Vector<T, dim> alpha_dd = coeff[0] * left_normal + coeff[1] * right_normal;
					space_der.col(1) = alpha_dd;
					Eigen::Vector<T, dim> L = std::pow(param_der(0, 0), 2) * left_ders(2, 0) + 2 * param_der(0, 0) * param_der(1, 0) * left_ders(1, 1) + std::pow(param_der(1, 0), 2) * left_ders(0, 2);

					Eigen::Vector<T, dim> temp = alpha_dd - L;
					dot_ders[0] = temp.dot(left_ders(1, 0));
					dot_ders[1] = temp.dot(left_ders(0, 1));
					param_der.template block<2, 1>(0, 1) = matSvd.solve(dot_ders);
					if (matSvd.info() != Eigen::Success)
						return ENUM_NURBS::NURBS_ERROR;

					L = std::pow(param_der(2, 0), 2) * right_ders(2, 0) + 2 * param_der(2, 0) * param_der(3, 0) * right_ders(1, 1) + std::pow(param_der(3, 0), 2) * right_ders(0, 2);
					temp = alpha_dd - L;
					dot_ders[0] = temp.dot(right_ders(1, 0));
					dot_ders[1] = temp.dot(right_ders(0, 1));
					param_der.template block<2, 1>(2, 1) = matSvd2.solve(dot_ders);
					if (matSvd2.info() != Eigen::Success)
						return ENUM_NURBS::NURBS_ERROR;
				}
				param_ders.push_back(param_der);
			}

			return ENUM_NURBS::NURBS_SUCCESS;
		}

		ENUM_NURBS estimate_next_param(const Eigen::Vector<T, 4>& initial_param, const Box<T, 4> domain, T& step, Eigen::Vector<T, 4>& next_param, int type)
		{
			int loop_count = 8;
			while (loop_count-- > 0)
			{
				//domain = Box<T, 4>(Eigen::Vector<T, 4>(0, 0, 0, 0), Eigen::Vector<T, 4>(1, 1, 1, 1));
				//Eigen::Vector<T, dim> tangent;
				Eigen::Vector<T, 4> k1;
				//Eigen::Vector<T, 4> k11;
				//Eigen::Vector<T, 4> k21;
				//Eigen::Vector<T, 4> k31;
				//Eigen::Vector<T, 4> k41;

				std::vector<Eigen::Matrix<T, 4, 1>> param_ders;
				std::vector<Eigen::Matrix<T, dim, 1>> space_ders;
				eval_preiamge_and_space_ders2(initial_param, param_ders, space_ders, type);
				// eval_preiamge_and_space_ders<1>(initial_param, param_ders, space_ders, type);
				k1 = step * param_ders[0];
				//eval_tangent(initial_param, tangent, k11, type);
				//k11 -= param_ders[0];
				//T dis1 = k11.norm();
				//k1 *= step;

				Eigen::Vector<T, 4> mid_param = initial_param + 0.5 * k1;
				if (domain.is_contain_point(mid_param) == false)
				{
					step /= 1.2;
					continue;
					//return ENUM_NURBS::NURBS_ERROR;
				}
				Eigen::Vector<T, 4> k2;
				eval_preiamge_and_space_ders2(mid_param, param_ders, space_ders, type);
				// eval_preiamge_and_space_ders<1>(mid_param, param_ders, space_ders, type);
				k2 = step * param_ders[0];
				//eval_tangent(mid_param, tangent, k21, type);
				//k21 -= param_ders[0];
				//T dis2 = k21.norm();
				//k2 *= step;
				mid_param = initial_param + 0.5 * k2;
				if (domain.is_contain_point(mid_param) == false)
				{
					step /= 1.2;
					continue;
					//return ENUM_NURBS::NURBS_ERROR;
				}
				Eigen::Vector<T, 4> k3;
				eval_preiamge_and_space_ders2(mid_param, param_ders, space_ders, type);
				// eval_preiamge_and_space_ders<1>(mid_param, param_ders, space_ders, type);
				k3 = step * param_ders[0];
				//eval_tangent(mid_param, tangent, k31, type);
				//k31 -= param_ders[0];
				//T dis3 = k31.norm();
				//k3 *= step;
				mid_param = initial_param + k3;
				if (domain.is_contain_point(mid_param) == false)
				{
					step /= 1.2;
					continue;
					//return ENUM_NURBS::NURBS_ERROR;
				}
				Eigen::Vector<T, 4> k4;
				eval_preiamge_and_space_ders2(mid_param, param_ders, space_ders, type);
				// eval_preiamge_and_space_ders<1>(mid_param, param_ders, space_ders, type);
				k4 = step * param_ders[0];
				//eval_tangent(mid_param, tangent, k41, type);
				//k41 -= param_ders[0];
				//T dis4 = k41.norm();
				//k4 *= step;

				next_param = initial_param + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
				if (domain.is_contain_point(next_param) == false)
				{
					step /= 1.2;
					continue;
					//return ENUM_NURBS::NURBS_ERROR;
				}
				return ENUM_NURBS::NURBS_SUCCESS;
			}
			return ENUM_NURBS::NURBS_ERROR;
		}


		bool may_contain_param(const Eigen::Vector2<T>& current_param, const Eigen::Vector2<T>& next_param, const Eigen::Vector2<T>& test_param)
		{
			Eigen::Vector2<T> param_vec = next_param - current_param;
			T vec_len = param_vec.norm();
			param_vec /= vec_len;
			//param_vec.normalize();
			Eigen::Vector2<T> test_vec = test_param - current_param;
			T project_len = test_vec.dot(param_vec);
			T eps = std::min(PRECISION<T>::value, TDEFAULT_ERROR<T>::value * vec_len);
			if (project_len > -eps && project_len < vec_len + eps)
			{
				return true;
			}
			return false;
		}


		bool is_point_on_arc(Eigen::Vector4<T>& current_param, Eigen::Vector2<T>& next_param, const Eigen::Vector2<T>& test_param, const Box<T, 4>& priori_enclosure, int type, bool is_positive)
		{
			Eigen::Vector<T, dim> test_point;
			m_right_surf_compute->point_on_surface(m_right_ders(0, 0), test_param[0], test_param[1], test_point);
			// m_right_ders(0, 0)->point_on_surface(test_param[0], test_param[1], test_point);
			Eigen::Vector<T, dim> start_point, end_point;
			m_right_surf_compute->point_on_surface(m_right_ders(0, 0), current_param[2], current_param[3], start_point);
			m_right_surf_compute->point_on_surface(m_right_ders(0, 0), next_param[0], next_param[1], end_point);
			// m_right_ders(0, 0)->point_on_surface(current_param[2], current_param[3], start_point);
			// m_right_ders(0, 0)->point_on_surface(next_param[0], next_param[1], end_point);
			T tolerance = type == 0 ? 1e-8 : 1e-4;
			if ((test_point - start_point).norm() < tolerance || (test_point - end_point).norm() < tolerance)
			{
				return true;
			}
			if (may_contain_param(current_param.template block<2, 1>(2, 0), next_param, test_param))
			{
				for (int iter_index = 0; iter_index < 2 * SURFACE_ITERATE_DEEP; ++iter_index)
				{

					Eigen::Vector<T, dim> chord_vec = end_point - start_point;
					T chord_length = chord_vec.norm();
					chord_vec /= chord_length;
					Eigen::Vector<T, dim> test_vec = test_point - start_point;
					T test_length = test_vec.dot(chord_vec);
					if (is_positive == false)
					{
						test_length = -test_length;
					}
					Eigen::Vector4<T> next_init_param, intersect_param;
					ENUM_NURBS error_code = estimate_next_param(current_param, priori_enclosure, test_length, next_init_param, type);
					if (error_code != ENUM_NURBS::NURBS_SUCCESS)
					{
						std::cout << "error " << std::endl;
						return false;
					}
					intersect_point_iteration(m_product_box, next_init_param, intersect_param);
					Eigen::Vector<T, dim> intersec_point;

					m_right_surf_compute->point_on_surface(m_right_ders(0, 0), intersect_param[2], intersect_param[3], intersec_point);
					// m_right_ders(0, 0)->point_on_surface(intersect_param[2], intersect_param[3], intersec_point);
					//double dis = (intersec_point - test_point).norm();
					if (priori_enclosure.is_contain_point(intersect_param) && (intersec_point - test_point).norm() < tolerance)
					{
						return true;
					}
					Eigen::Vector2<T> mid_param = intersect_param.template block<2, 1>(2, 0);

					bool flag1 = may_contain_param(current_param.template block<2, 1>(2, 0), mid_param, test_param);
					bool flag2 = may_contain_param(mid_param, next_param, test_param);

					if (flag1 == true && flag2 == true)
					{
						std::cout << "error" << std::endl;
					}
					if (flag1 == true)
					{
						next_param = mid_param;
						m_right_surf_compute->point_on_surface(m_right_ders(0, 0), next_param[0], next_param[1], end_point);
						// m_right_ders(0, 0)->point_on_surface(next_param[0], next_param[1], end_point);
					}
					else if (flag2 == true)
					{
						current_param.template block<2, 1>(2, 0) = mid_param;
						m_right_surf_compute->point_on_surface(m_right_ders(0, 0), current_param[2], current_param[3], start_point);
						// m_left_surf_compute->find_point_on_surface(intersec_point, current_param[0], current_param[1]);
						// m_right_ders(0, 0)->point_on_surface(current_param[2], current_param[3], start_point);
						m_left_ders(0, 0)->find_point_on_surface(intersec_point, current_param[0], current_param[1]);
					}
					else
					{
						break;
					}

				}
			}
			return false;
		}

		bool is_closed(surfs_int_points_chat<T, dim>& int_chat, const surf_surf_intersect_point<T, dim>& int_point, int type)
		{
			// TODO: 重写
			if (int_chat.m_inter_points.size() > 3 && int_point.m_priori_enclosure.is_contain_point(int_chat.m_inter_points[0].m_uv))
			{
				Eigen::Vector2<T> test_param = int_chat.m_inter_points[0].m_uv.template block<2, 1>(2, 0);
				Eigen::Vector2<T> next_param = int_point.m_uv.template block<2, 1>(2, 0);
				Eigen::Vector4<T> current_param = int_chat.m_inter_points.back().m_uv;
				return is_point_on_arc(current_param, next_param, test_param, int_point.m_priori_enclosure, type, int_chat.m_is_positive_direction);
			}
			return false;
		}

		bool is_point_on_intcurve(const Eigen::Vector4<T>& param, bool is_transversal, int type)
		{
			Eigen::Vector2<T> test_param = param.template block<2, 1>(2, 0);
			for (surfs_int_points_chat<T, dim>& int_chat : m_result.m_int_chats)
			{
				if (is_transversal != int_chat.m_is_transversal)
				{
					continue;
				}
				int int_points_count = int_chat.m_inter_points.size();
				for (int index = 0; index < int_points_count - 1; ++index)
				{
					const Box<T, 4>& priori_enclosure = int_chat.m_inter_points[index].m_priori_enclosure;
					if (priori_enclosure.is_contain_point(param))
					{
						Eigen::Vector4<T> current_param = int_chat.m_inter_points[index].m_uv;
						//Eigen::Vector2<T> next_param = int_chat.m_inter_points[index + 1].m_uv.template block<2, 1>(0, 0);
						Eigen::Vector2<T> next_param = int_chat.m_inter_points[index + 1].m_uv.template block<2, 1>(2, 0);

						//Eigen::Vector4<T> next_param = int_chat.m_inter_points[index + 1].m_uv;
						if (is_point_on_arc(current_param, next_param, test_param, priori_enclosure, type, int_chat.m_is_positive_direction) == true)
						{
							return true;
						}

					}
				}
			}
			return false;
		}


		ENUM_NURBS eval_priori_enclosure_inner(const Box<T, 4>& initial_box, const Box<T, 4>& param_box, T step_size, Box<T, 4>& priori_enclosure, int& type, bool& may_arrived_singular, bool transversal = true)
		{
			may_arrived_singular = false;
			std::array<Box<T, 4>, 1> param_ders;
			std::array<Box<T, dim>, 1> space_ders;
			ENUM_NURBS error_code;
			if (transversal == true)
			{
				error_code = eval_preiamge_and_space_ders<1>(param_box, param_ders, space_ders, may_arrived_singular);
				type = 0;
			}
			else
			{
				error_code = eval_preiamge_and_space_ders_tangent<1>(param_box, param_ders, space_ders, type);
			}
			if (ENUM_NURBS::NURBS_SUCCESS != error_code)
			{
				return error_code;
			}

			Interval<T> step_box(0, step_size);

			Box<T, 4> box_interval = param_ders[0] * step_box;
			Box<T, 4> temp_priori_enclosure = initial_box + box_interval;

			if (true == m_product_box.is_contain_box(temp_priori_enclosure, 5.0 * PRECISION<T>::value))
			{
				temp_priori_enclosure.intersect(m_product_box, priori_enclosure);
				return ENUM_NURBS::NURBS_SUCCESS;
			}
			return ENUM_NURBS::NURBS_ERROR;


		}

		ENUM_NURBS eval_priori_enclosure_inner(const Eigen::Vector4<T>& initial_param, const Box<T, 4>& param_box, T step_size, Box<T, 4>& priori_enclosure, int& type, bool& may_arrived_singular, bool transversal = true)
		{
			may_arrived_singular = false;
			std::array<Box<T, 4>, 1> param_ders;
			std::array<Box<T, dim>, 1> space_ders;
			ENUM_NURBS error_code;
			if (transversal == true)
			{
				error_code = eval_preiamge_and_space_ders<1>(param_box, param_ders, space_ders, may_arrived_singular);
				type = 0;
			}
			else
			{
				error_code = eval_preiamge_and_space_ders_tangent<1>(param_box, param_ders, space_ders, type);
			}
			if (ENUM_NURBS::NURBS_SUCCESS != error_code)
			{
				return error_code;
			}

			Interval<T> step_box(0, step_size);

			Box<T, 4> temp_priori_enclosure = param_ders[0] * step_box;
			temp_priori_enclosure.move_vector(initial_param);
			// Box<T, 4> temp_priori_enclosure = initial_box + box_interval;

			if (true == m_product_box.is_contain_box(temp_priori_enclosure, 5.0 * PRECISION<T>::value))
			{
				temp_priori_enclosure.intersect(m_product_box, priori_enclosure);
				return ENUM_NURBS::NURBS_SUCCESS;
			}
			return ENUM_NURBS::NURBS_ERROR;


		}
		
		bool eval_priori_enclosure(const Box<T, 4>& initial_box, const T& bigger_step_size, T& smalll_step_size, Box<T, 4>& priori_enclosure, T& angle_diff, int& type, bool& may_arrived_singular, bool transversal = true)
		{
			may_arrived_singular = false;
			int loopCount = 2;
			angle_diff = 4.0;
			if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, initial_box, bigger_step_size, priori_enclosure, type, may_arrived_singular, transversal))
			{
				angle_diff = 4.0;
				return false;
			}

			while (loopCount > 0)
			{
				loopCount -= 1;
				//if (loopCount == 4)
				//{
				//    Eigen::Vector<T, 4> &min = priori_enclosure.Min;
				//    Eigen::Vector<T, 4> &max = priori_enclosure.Max;

				//    for (int index = 0; index < 2; ++index)
				//    {
				//        T startIndex = 2 * index;
				//        if ((max[startIndex] - min[startIndex]) > 10 * (max[startIndex + 1] - min[startIndex + 1]))
				//        {
				//            T mid = (max[1 + startIndex] + min[1 + startIndex]) / 2.0;
				//            T len = (max[startIndex] - min[startIndex]) / 10;
				//            max[1 + startIndex] = mid + len;
				//            min[1 + startIndex] = mid - len;
				//        }
				//        else if ((max[1 + startIndex] - min[1 + startIndex]) > 10 * (max[startIndex] - min[startIndex]))
				//        {
				//            T mid = (max[startIndex] + min[startIndex]) / 2.0;
				//            T len = (max[1 + startIndex] - min[1 + startIndex]) / 10;
				//            max[startIndex] = mid + len;
				//            min[startIndex] = mid - len;
				//        }
				//    }
				//    priori_enclosure.intersect(m_product_box, priori_enclosure);
				//}
				Box<T, 4> next_priori_enclosure;
				smalll_step_size /= 2.0;
				if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, priori_enclosure, smalll_step_size, next_priori_enclosure, type, may_arrived_singular, transversal))
				{
					// Box<T, 2> left_param_box(priori_enclosure.Min.template block<2, 1>(0, 0), priori_enclosure.Max.template block<2, 1>(0, 0));
					// Box<T, 2> right_param_box(priori_enclosure.Min.template block<2, 1>(2, 0), priori_enclosure.Max.template block<2, 1>(2, 0));
					// Box<T, dim> left_normal_box, right_normal_box;
					// if constexpr (left_surface_type::is_ratio == false)
					// {
					// 	//使用精确的normal_surface和使用两个box叉乘得到的结果竟然基本一样？？？？？？？
					// 	m_left_normal_surface->get_box(left_param_box, left_normal_box);
					// }
					// else
					// {
					// 	Eigen::MatrixX<Box<T, dim>> left_boxes;

					// 	ENUM_NURBS error_code = bounding_box<left_surface_type>::get_ders_bounding_box(1, m_left_ders, left_param_box, left_boxes);
					// 	if (error_code != ENUM_NURBS::NURBS_SUCCESS)
					// 	{
					// 		angle_diff = M_PI;
					// 		return false;
					// 	}
					// 	left_normal_box = interval_algorithm::cross(left_boxes(1, 0), left_boxes(0, 1));
					// }
					// if constexpr (right_surface_type::is_ratio == false)
					// {
					// 	//使用精确的normal_surface和使用两个box叉乘得到的结果竟然基本一样？？？？？？？
					// 	m_right_normal_surface->get_box(right_param_box, right_normal_box);
					// }
					// else
					// {
					// 	Eigen::MatrixX<Box<T, dim>> right_boxes;

					// 	ENUM_NURBS error_code = bounding_box<right_surface_type>::get_ders_bounding_box(1, m_right_ders, right_param_box, right_boxes);
					// 	if (error_code != ENUM_NURBS::NURBS_SUCCESS)
					// 	{
					// 		angle_diff = M_PI;
					// 		return false;
					// 	}
					// 	right_normal_box = interval_algorithm::cross(right_boxes(1, 0), right_boxes(0, 1));
					// }

					// Eigen::Vector<T, dim> origin;
					// origin.setConstant(0.0);
					// cone<T, dim> normal_cone = point_box(left_normal_box, origin);
					// angle_diff = normal_cone.m_angle;
					// normal_cone = point_box(right_normal_box, origin);
					// angle_diff = std::max(angle_diff, normal_cone.m_angle);

					return false;
				}
				if (priori_enclosure.is_contain_box(next_priori_enclosure))
				{
					//if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, next_priori_enclosure, smalll_step_size, priori_enclosure, type, may_arrived_singular, transversal))
					//{
					//	return false;
					//}
					return true;
				}
				
				//else
				if (loopCount == 1)
				{
					Box<T, 4> temp_box = next_priori_enclosure;
					if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, temp_box, smalll_step_size, next_priori_enclosure, type, may_arrived_singular, transversal))
					{
						// priori_enclosure = temp_box;
						// Box<T, 2> left_param_box(priori_enclosure.Min.template block<2, 1>(0, 0), priori_enclosure.Max.template block<2, 1>(0, 0));
						// Box<T, 2> right_param_box(priori_enclosure.Min.template block<2, 1>(2, 0), priori_enclosure.Max.template block<2, 1>(2, 0));
						// Box<T, dim> left_normal_box, right_normal_box;
						// if constexpr (left_surface_type::is_ratio == false)
						// {
						// 	//使用精确的normal_surface和使用两个box叉乘得到的结果竟然基本一样？？？？？？？
						// 	m_left_normal_surface->get_box(left_param_box, left_normal_box);
						// }
						// else
						// {
						// 	Eigen::MatrixX<Box<T, dim>> left_boxes;
						// 	ENUM_NURBS error_code = bounding_box<left_surface_type>::get_ders_bounding_box(1, m_left_ders, left_param_box, left_boxes);
						// 	if (error_code != ENUM_NURBS::NURBS_SUCCESS)
						// 	{
						// 		angle_diff = M_PI;
						// 		return false;
						// 	}

						// 	left_normal_box = interval_algorithm::cross(left_boxes(1, 0), left_boxes(0, 1));
						// }

						// if constexpr (right_surface_type::is_ratio == false)
						// {
						// 	//使用精确的normal_surface和使用两个box叉乘得到的结果竟然基本一样？？？？？？？
						// 	m_right_normal_surface->get_box(right_param_box, right_normal_box);
						// }
						// else
						// {
						// 	Eigen::MatrixX<Box<T, dim>> right_boxes;
						// 	ENUM_NURBS error_code = bounding_box<right_surface_type>::get_ders_bounding_box(1, m_right_ders, right_param_box, right_boxes);
						// 	if (error_code != ENUM_NURBS::NURBS_SUCCESS)
						// 	{
						// 		angle_diff = M_PI;
						// 		return false;
						// 	}

						// 	right_normal_box = interval_algorithm::cross(right_boxes(1, 0), right_boxes(0, 1));
						// }
						// Eigen::Vector<T, dim> origin;
						// origin.setConstant(0.0);
						// cone<T, dim> normal_cone = point_box(left_normal_box, origin);
						// angle_diff = normal_cone.m_angle;
						// normal_cone = point_box(right_normal_box, origin);
						// angle_diff = std::max(angle_diff, normal_cone.m_angle);

						return false;
					}
					if (temp_box.is_contain_box(next_priori_enclosure))
					{
						priori_enclosure = temp_box;
						//if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_box, next_priori_enclosure, smalll_step_size, priori_enclosure, type, may_arrived_singular, transversal))
						//{
						//	return false;
						//}
						return true;
					}
					priori_enclosure = next_priori_enclosure;

				}
				//smalll_step_size /= 1.5;
				if (std::abs(smalll_step_size) < TDEFAULT_ERROR<T>::value)
				{
					return false;
				}
			}
			return false;
		}

		
		bool eval_priori_enclosure(const Eigen::Vector4<T>& initial_param, const Eigen::Vector4<T>& param_tangent_vector, const T& bigger_step_size, T& smalll_step_size, Box<T, 4>& priori_enclosure, T& angle_diff, int& type, bool& may_arrived_singular, bool transversal = true)
		{
			may_arrived_singular = false;
			int loopCount = 2;
		
			T sacle = 1.0;

			Eigen::Vector4<T> end_param = initial_param + bigger_step_size * param_tangent_vector;
			end_param = initial_param + (bigger_step_size * sacle) * param_tangent_vector;
			for (int i = 0; i < 4; ++i)
			{
				priori_enclosure.Min[i] = std::min(initial_param[i], end_param[i]);
				priori_enclosure.Max[i] = std::max(initial_param[i], end_param[i]);
				priori_enclosure.Max[i] += 1e-6;
				priori_enclosure.Min[i] -= 1e-6;
			}

			// TODO: 
			priori_enclosure.intersect(m_product_box, priori_enclosure);
			if ((priori_enclosure.Max - priori_enclosure.Min).norm() < 1e-8) 
			{
				smalll_step_size /= 2.0;
				return false;
			}

			while (loopCount > 0)
			{
				loopCount -= 1;
				Box<T, 4> next_priori_enclosure;
				smalll_step_size /= 2.0;
				if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_param, priori_enclosure, smalll_step_size, next_priori_enclosure, type, may_arrived_singular, transversal))
				{
					return false;
				}
				if (priori_enclosure.is_contain_box(next_priori_enclosure))
				{
					return true;
				}
				
				//else
				if (loopCount == 1)
				{
					Box<T, 4> temp_box = next_priori_enclosure;
					if (ENUM_NURBS::NURBS_SUCCESS != eval_priori_enclosure_inner(initial_param, temp_box, smalll_step_size, next_priori_enclosure, type, may_arrived_singular, transversal))
					{
						return false;
					}
					if (temp_box.is_contain_box(next_priori_enclosure))
					{
						priori_enclosure = temp_box;
						return true;
					}
					priori_enclosure = next_priori_enclosure;

				}
				//smalll_step_size /= 1.5;
				if (std::abs(smalll_step_size) < TDEFAULT_ERROR<T>::value)
				{
					return false;
				}
			}
			return false;
		}

		//迭代加细交点(周期性曲面待处理)
		ENUM_NURBS intersect_point_iteration(const Box<T, 4>& domian, Eigen::Vector<T, 4> current_param, Eigen::Vector<T, 4>& intersect_param)
		{
			// bool left_surf_u_closed = left_surf->is_u_closed();
			// bool left_surf_v_closed = left_surf->is_v_closed();
			// bool right_surf_u_closed = right_surf->is_u_closed();
			// bool right_surf_v_closed = right_surf->is_v_closed();

			// Eigen::Vector<T, dim> left_point, right_point;
			Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> left_ders, right_ders;

			Eigen::Vector<T, dim> vec;
			T min_distance = 100000;
			intersect_param = current_param;

			//迭代次数需要修改
			for (int loop_index = 0; loop_index < 2 * SURFACE_ITERATE_DEEP; ++loop_index)
			{
				m_left_surf_compute->template derivative_on_surface<1>(current_param[0], current_param[1], left_ders);
				m_right_surf_compute->template derivative_on_surface<1>(current_param[2], current_param[3], right_ders);
				// m_left_ders(0, 0)->template derivative_on_surface<1>(current_param[0], current_param[1], left_ders);
				// m_right_ders(0, 0)->template derivative_on_surface<1>(current_param[2], current_param[3], right_ders);
				// m_left_surface->point_on_surface(current_param[0], current_param[1], left_point);
				// m_right_surface->point_on_surface(current_param[2], current_param[3], right_point);
				// vec = right_point - left_point;
				vec = right_ders(0, 0) - left_ders(0, 0);
				T distance = vec.squaredNorm();
				if (distance < min_distance)
				{
					min_distance = distance;
					intersect_param = current_param;

					bool flag = true;
					for (int index = 0; index < dim; ++index)
					{
						if (std::abs(vec[index]) > PRECISION<T>::value)
						{
							flag = false;
							break;
						}
					}
					if (flag == true)
						return ENUM_NURBS::NURBS_SUCCESS;
				}

				Eigen::Matrix<T, dim, 4> mat;
				// Eigen::Vector<T, dim> temp;
				// m_left_u_tangent_surface->point_on_surface(current_param[0], current_param[1], temp);
				// mat.col(0) = temp;
				mat.col(0) = left_ders(1, 0);
				// m_left_v_tangent_surface->point_on_surface(current_param[0], current_param[1], temp);
				// mat.col(1) = temp;
				mat.col(1) = left_ders(0, 1);
				// m_right_u_tangent_surface->point_on_surface(current_param[2], current_param[3], temp);
				// mat.col(2) = -1.0 * temp;
				mat.col(2) = -1.0 * right_ders(1, 0);
				// m_right_v_tangent_surface->point_on_surface(current_param[2], current_param[3], temp);
				// mat.col(3) = -1.0 * temp;
				mat.col(3) = -1.0 * right_ders(0, 1);
				Eigen::JacobiSVD<Eigen::Matrix<T, dim, 4>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
				Eigen::Vector<T, 4> delta = matSvd.solve(vec);
				if (matSvd.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;

				Eigen::Vector<T, 4> next_param = current_param + delta;
				for (int index = 0; index < 4; ++index)
				{
					if (next_param[index] < domian.Min[index])
						next_param[index] = domian.Min[index];
					else if (next_param[index] > domian.Max[index])
						next_param[index] = domian.Max[index];
				}


				// bool is_closed_flag = cur->is_closed();
				// if (left_surf_u_closed)
				// {
				//     if (next_u < min)
				//         next_u = max - (min - next_u);
				//     else if (next_u > max)
				//         next_u = min + (next_u - max);
				// }
				// else
				// {
				//     if (next_u < min)
				//         next_u = min;
				//     else if (next_u > max)
				//         next_u = max;
				// }

				current_param = next_param;



			}

			return ENUM_NURBS::NURBS_ERROR;
		}

		//迭代加细交点(周期性曲面待处理)
		ENUM_NURBS intersect_singular_point_iteration(const Box<T, 4>& domian, Eigen::Vector<T, 4> current_param, Eigen::Vector<T, 4>& intersect_param)
		{
			intersect_param = current_param;
			T min_distance = 10000000;
			//迭代次数需要修改
			for (int loop_index = 0; loop_index < SURFACE_ITERATE_DEEP; ++loop_index)
			{
				Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_derivatives, right_derivatives;
				m_left_surf_compute->template derivative_on_surface<2>(current_param[0], current_param[1], left_derivatives);
				m_right_surf_compute->template derivative_on_surface<2>(current_param[2], current_param[3], right_derivatives);
				// m_left_ders(0, 0)->template derivative_on_surface<2>(current_param[0], current_param[1], left_derivatives);
				// m_right_ders(0, 0)->template derivative_on_surface<2>(current_param[2], current_param[3], right_derivatives);
				Eigen::Vector<T, dim + 1> vec;
				vec.template block<dim, 1>(0, 0) = right_derivatives(0, 0) - left_derivatives(0, 0);


				Eigen::Matrix<T, dim + 1, 4> mat;
				mat.template block<dim, 1>(0, 0) = left_derivatives(1, 0);
				mat.template block<dim, 1>(0, 1) = left_derivatives(0, 1);
				mat.template block<dim, 1>(0, 2) = -right_derivatives(1, 0);
				mat.template block<dim, 1>(0, 3) = -right_derivatives(0, 1);

				Eigen::Vector<T, dim> dsdt = right_derivatives(1, 0).cross(right_derivatives(0, 1));
				Eigen::Vector<T, dim> dudv = left_derivatives(1, 0).cross(left_derivatives(0, 1));
				Eigen::Vector<T, dim> dudvdsdt = dudv.cross(dsdt);
				Eigen::Vector<T, dim> dtdss = right_derivatives(0, 1).cross(right_derivatives(2, 0));
				Eigen::Vector<T, dim> dsdst = right_derivatives(1, 0).cross(right_derivatives(1, 1));
				Eigen::Vector<T, dim> dtdst = right_derivatives(0, 1).cross(right_derivatives(1, 1));
				Eigen::Vector<T, dim> dsdtt = right_derivatives(1, 0).cross(right_derivatives(0, 2));

				Eigen::Vector<T, dim> dvduu = left_derivatives(0, 1).cross(left_derivatives(2, 0));
				Eigen::Vector<T, dim> duduv = left_derivatives(1, 0).cross(left_derivatives(1, 1));
				Eigen::Vector<T, dim> dvduv = left_derivatives(0, 1).cross(left_derivatives(1, 1));
				Eigen::Vector<T, dim> dudvv = left_derivatives(1, 0).cross(left_derivatives(0, 2));

				double coeff = 2.0 / (dudv.squaredNorm() * dsdt.squaredNorm());
				//double coeff = 1.0;
				mat(dim, 0) = coeff * (-dvduu + duduv).cross(dsdt).dot(dudvdsdt);
				mat(dim, 1) = coeff * (-dvduv + dudvv).cross(dsdt).dot(dudvdsdt);
				mat(dim, 2) = coeff * dudv.cross(-dtdss + dsdst).dot(dudvdsdt);
				mat(dim, 3) = coeff * dudv.cross(-dtdst + dsdtt).dot(dudvdsdt);

				vec[dim] = -coeff * dudvdsdt.squaredNorm();
				T distance = vec.squaredNorm();
				if (distance < min_distance)
				{
					min_distance = distance;
					intersect_param = current_param;

					bool flag = true;
					for (int index = 0; index <= dim; ++index)
					{
						if (std::abs(vec[index]) > PRECISION<T>::value)
						{
							flag = false;
							break;
						}
					}
					if (flag == true)
						return ENUM_NURBS::NURBS_SUCCESS;
				}

				Eigen::JacobiSVD<Eigen::Matrix<T, dim + 1, 4>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
				Eigen::Vector<T, 4> delta = matSvd.solve(vec);
				if (matSvd.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;

				Eigen::Vector<T, 4> next_param = current_param + delta;
				for (int index = 0; index < 4; ++index)
				{
					if (next_param[index] < domian.Min[index])
						next_param[index] = domian.Min[index];
					else if (next_param[index] > domian.Max[index])
						next_param[index] = domian.Max[index];
				}

				current_param = next_param;

			}

			return ENUM_NURBS::NURBS_ERROR;
		}

		ENUM_NURBS eval_singular_points(const Eigen::Vector<T, 4>& intersect_param, const std::vector<Eigen::Matrix<T, 4, 1>>& param_ders, const std::vector<Eigen::Matrix<T, dim, 1>>& space_ders,
			std::vector<surfs_int_points_chat<T, dim>>& singular_chats)
		{
			for (const Box<T, 4>&singular_box : m_singular_boxes)
			{
				if (singular_box.is_contain_point(intersect_param, TDEFAULT_ERROR<T>::value * 0.01))
				{
					return ENUM_NURBS::NURBS_SUCCESS;
				}
			}
			Box<T, 4> singular_box;
			create_box(intersect_param, TDEFAULT_ERROR<T>::value * 0.01, singular_box);
			m_singular_boxes.push_back(singular_box);
			surf_surf_intersect_point<T, dim> singluar_point;
			singluar_point.m_uv = intersect_param;
			m_left_surf_compute->point_on_surface(m_left_ders(0, 0), intersect_param[0], intersect_param[1], singluar_point.m_point);
			// m_left_ders(0, 0)->point_on_surface(intersect_param[0], intersect_param[1], singluar_point.m_point);
			Eigen::Vector<T, 4> translate(TDEFAULT_ERROR<T>::value, TDEFAULT_ERROR<T>::value, TDEFAULT_ERROR<T>::value, TDEFAULT_ERROR<T>::value);
			singluar_point.m_priori_enclosure = Box<T, 4>(intersect_param - translate, intersect_param + translate);
			std::vector<Eigen::Vector4<T>> next_params(4);
			for (int index = 0; index < 2; ++index)
			{
				Eigen::Vector4<T> param_tangent = param_ders[index].col(0).normalized();
				Eigen::Vector4<T> initial_param1 = intersect_param + param_tangent * TDEFAULT_ERROR<T>::value;
				Eigen::Vector4<T> intial_param2 = intersect_param - param_tangent * TDEFAULT_ERROR<T>::value;
				intersect_point_iteration(m_product_box, initial_param1, next_params[2 * index]);
				intersect_point_iteration(m_product_box, intial_param2, next_params[2 * index + 1]);
			}

			Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_ders;
			// m_left_ders(0, 0)->template derivative_on_surface<2, 2>(intersect_param[0], intersect_param[1], left_ders);
			m_left_surf_compute->template derivative_on_surface<2>(intersect_param[0], intersect_param[1], left_ders);
			Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
			left_normal.normalize();
			Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_ders;
			// m_right_ders(0, 0)->template derivative_on_surface<2, 2>(intersect_param[2], intersect_param[3], right_ders);
			m_right_surf_compute->template derivative_on_surface<2>(intersect_param[2], intersect_param[3], right_ders);
			Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
			right_normal.normalize();
			bool flag = left_normal.dot(right_normal) > 0;

			std::vector<Eigen::Matrix<T, 4, 1>> current_param_ders;
			std::vector<Eigen::Matrix<T, dim, 1>> current_space_ders;
			for (int index = 0; index < 4; ++index)
			{
				singular_chats.push_back(surfs_int_points_chat<T, dim>());
				surfs_int_points_chat<T, dim>& end_chat = singular_chats.back();
				std::vector<surf_surf_intersect_point<T, dim>>& int_point = end_chat.m_inter_points;
				int_point.resize(2);
				int_point[0] = singluar_point;

				Eigen::Vector<T, 4> singular_intersect_param;
				std::array<T, 4> ks;
				std::array<T, 2> left_EG, right_EG;
				T k1, k2, k3, k4;
				if (intersect_singular_point_iteration(m_product_box, end_chat.m_inter_points[1].m_uv, singular_intersect_param) == ENUM_NURBS::NURBS_SUCCESS)
				{
					eval_preiamge_and_space_ders<1>(singular_intersect_param, current_param_ders, current_space_ders, ks, left_EG, right_EG, 3);
					if (current_space_ders.size() == 1)
					{
						end_chat.m_is_transversal = false;
						end_chat.m_inter_points[1].m_uv = singular_intersect_param;
						// m_left_ders(0, 0)->point_on_surface(singular_intersect_param[0], singular_intersect_param[1], end_chat.m_inter_points[1].m_point);
						m_left_surf_compute->point_on_surface(m_left_ders(0, 0), singular_intersect_param[0], singular_intersect_param[1], end_chat.m_inter_points[1].m_point);
						end_chat.start_point_type = POINT_TYPE::SINGULAR;

						T tangent_dot = current_param_ders[0].dot(param_ders[index / 2]);
						if (index % 2 == 1)
						{
							tangent_dot *= -1;
						}
						end_chat.m_is_positive_direction = tangent_dot > 0;
						continue;
					}
				}

				end_chat.m_is_transversal = true;
				end_chat.start_point_type = POINT_TYPE::SINGULAR;
				int_point[1].m_uv = next_params[index];
				// m_left_ders(0, 0)->point_on_surface(next_params[index][0], next_params[index][1], int_point[1].m_point);
				m_left_surf_compute->point_on_surface(m_left_ders(0, 0), next_params[index][0], next_params[index][1], int_point[1].m_point);

				eval_preiamge_and_space_ders<1>(end_chat.m_inter_points[1].m_uv, current_param_ders, current_space_ders, ks, left_EG, right_EG, 0);
				T tangent_dot = current_param_ders[0].dot(param_ders[index / 2]);
				if (index % 2 == 1)
				{
					tangent_dot *= -1;
				}
				end_chat.m_is_positive_direction = tangent_dot > 0;
			}

			return ENUM_NURBS::NURBS_SUCCESS;
		}

		ENUM_NURBS try_find_singular_points(const Eigen::Vector<T, 4>& initial_param, const Box<T, 4>& priori_enclosure, std::vector<surfs_int_points_chat<T, dim>>& singular_chats)
		{
			// try find singular point
			singular_chats.clear();
			Eigen::Vector<T, 4> pre_param, intersect_param;
			if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_iteration(m_product_box, initial_param, pre_param))
			{
				if (ENUM_NURBS::NURBS_SUCCESS == intersect_singular_point_iteration(m_product_box, pre_param, intersect_param))
				{
					if (priori_enclosure.is_contain_point(intersect_param) == true)
					{
						// TODO: 整理
						std::vector<Eigen::Matrix<T, 4, 1>> param_ders_temp;
						std::vector<Eigen::Matrix<T, dim, 1>> space_ders_temp;
						ENUM_NURBS error_code = eval_preiamge_and_space_ders<1>(intersect_param, param_ders_temp, space_ders_temp, 3);
						if (space_ders_temp.size() == 2)
						{
							eval_singular_points(intersect_param, param_ders_temp, space_ders_temp, singular_chats);
						}
					}
				}
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		};


		T estimate_arc_length(const std::array<T, 4>& ks, const std::array<T, 2>& left_EG, const std::array<T, 2>& right_EG, const T curvature, const Eigen::Vector4<T>& param_tangent)
		{
			T arc_length = 2.0 * m_angle_eps;
			if (curvature < TDEFAULT_ERROR<T>::value)
			{
				arc_length = TINFINITE<T>::value;
			}
			else
			{
				arc_length /= curvature;
			}
			Eigen::Vector4<T> param_vector_len = param_tangent * arc_length;
			
			T x1 = std::abs(param_vector_len[0]);
			T x2 = std::abs(param_vector_len[1]);
			T x3 = std::abs(param_vector_len[2]);
			T x4 = std::abs(param_vector_len[3]);
			T scale = 1.0;
			T scale_angle = 1.0;
			if (std::abs(ks[0]) < TDEFAULT_ERROR<T>::value || x1 < TDEFAULT_ERROR<T>::value)
			{
				// arc_length = TINFINITE<T>::value;
			}
			else
			{
				scale = std::min(scale, std::abs(scale_angle * m_angle_eps / (ks[0] * std::sqrt(left_EG[0]))) / x1);
			}
			
			if (std::abs(ks[1]) < TDEFAULT_ERROR<T>::value || x2 < TDEFAULT_ERROR<T>::value)
			{
				// arc_length = TINFINITE<T>::value;
			}
			else
			{
				scale = std::min(scale, std::abs(scale_angle * m_angle_eps / (ks[1] * std::sqrt(left_EG[1]))) / x2);
			}

			if (std::abs(ks[2]) < TDEFAULT_ERROR<T>::value || x3 < TDEFAULT_ERROR<T>::value)
			{
				// arc_length = TINFINITE<T>::value;
			}
			else
			{
				scale = std::min(scale, std::abs(scale_angle * m_angle_eps / (ks[2] * std::sqrt(right_EG[0]))) / x3);
			}
			
			if (std::abs(ks[3]) < TDEFAULT_ERROR<T>::value || x4 < TDEFAULT_ERROR<T>::value)
			{
				// arc_length = TINFINITE<T>::value;
			}
			else
			{
				scale = std::min(scale, std::abs(scale_angle * m_angle_eps / (ks[3] * std::sqrt(right_EG[1]))) / x4);
			}

			arc_length *= scale;
			return arc_length;
		}

		ENUM_NURBS trace_point(surfs_int_points_chat<T, dim>& current_intpoints, bool& is_stop, /* bool& is_arrived_singular, */ T box_len = TDEFAULT_ERROR<T>::value * 0.01)
		{
			//bool direction;
			is_stop = false;
			bool transversal = current_intpoints.m_is_transversal;
			surf_surf_intersect_point<T, dim>& current_intpoint = current_intpoints.m_inter_points.back();
			Eigen::Vector<T, 4> param = current_intpoints.m_inter_points.back().m_uv;
			Box<T, 4> initial_box;
			create_box(param, box_len, initial_box);
			initial_box.intersect(m_product_box, initial_box);


			std::vector<Eigen::Matrix<T, 4, 2>> param_ders;
			std::vector<Eigen::Matrix<T, dim, 2>> space_ders;
			std::array<T, 4> ks;
			std::array<T, 2> left_EG, right_EG;
			eval_preiamge_and_space_ders<2>(param, param_ders, space_ders, ks, left_EG, right_EG, transversal == true ? 0 : 3);
			// TODO : 处理branch
			
			T curvature = space_ders[0].col(1).norm();
			Eigen::Vector4<T> param_tangent = param_ders[0].col(0);
			T arc_length = estimate_arc_length(ks, left_EG, right_EG, curvature, param_tangent);
			
			// Eigen::Vector4<T> param_tangent = param_ders[0].col(0);
			// TODO : 处理step太大超过参数域
			T bigger_step_size = current_intpoints.m_is_positive_direction == false ? -arc_length : arc_length;
			if (std::abs(bigger_step_size) > 10.0 * std::abs(current_intpoint.m_previous_step))
			{
				bigger_step_size = 10 * current_intpoint.m_previous_step;
			}
			T small_step = bigger_step_size;
			int type;
			T angle_diff;
			Box<T, 4> priori_enclosure;
			Eigen::Vector<T, 4> intersect_param;
			Eigen::Vector<T, 4> next_init_param;
			bool may_arrived_singular = false;
			int counttemp = looop_count;
			// std::cout << "loopCount_begin: " << looop_count << std::endl;
			// while (false == eval_priori_enclosure(initial_box, bigger_step_size, small_step, priori_enclosure, angle_diff, type, may_arrived_singular, transversal))
			// TODO: 处理
			T is_negetive = bigger_step_size > 0 ? 1 : -1;
			Eigen::Vector4<T> next_minial_point = param + is_negetive * 1e-4 * param_tangent;
			if (m_product_box.is_contain_point(next_minial_point, 1e-10) == false)
			{
				return ENUM_NURBS::NURBS_ERROR;
			}
			while (false == eval_priori_enclosure(param, param_tangent, bigger_step_size, small_step, priori_enclosure, angle_diff, type, may_arrived_singular, transversal))
			{
				// TODO: delete
				if (std::abs(bigger_step_size) < 1e-4)
				{
					return ENUM_NURBS::NURBS_ERROR;

				}

				// bigger_step_size /= 2.0;
				bigger_step_size = small_step;
				small_step = bigger_step_size;
			}
			
			// std::cout << "loopCount: " << looop_count << std::endl;
			// std::cout << "loopCountloop: " << looop_count - counttemp << std::endl;

			int inde_t = 0;
			do
			{
				if (ENUM_NURBS::NURBS_SUCCESS != estimate_next_param(param, priori_enclosure, small_step, next_init_param, type))
				{
					inde_t += 1;
					// Eigen::Vector3d ppp;
					// m_left_ders(0, 0)->point_on_surface(param[0], param[1], ppp);
					if (std::abs(small_step) < 1e-4)
					{
						return ENUM_NURBS::NURBS_ERROR;

					}
					continue;
				}
				if (ENUM_NURBS::NURBS_SUCCESS == intersect_point_iteration(m_product_box, next_init_param, intersect_param))
				{
					if (transversal == false)
					{
						next_init_param = intersect_param;
						if (ENUM_NURBS::NURBS_SUCCESS == intersect_singular_point_iteration(m_product_box, next_init_param, intersect_param))
						{
							if (priori_enclosure.is_contain_point(intersect_param) == true)
							{
								break;
							}
						}
					}
					else
					{
						if (priori_enclosure.is_contain_point(intersect_param) == true)
						{
							break;
						}
					}
				}
				small_step /= 1.2;

			} while (true);

			if (inde_t > 1)
			{
				std::cout << "inde_t: " << std::endl;
			}

			for (int index = 0; index < m_boundary_chats.size(); ++index)
			{
				surf_surf_intersect_point<T, dim>& boundary_point = m_boundary_chats[index].m_inter_points.back();

				if (priori_enclosure.is_contain_point(boundary_point.m_uv))
				{
					surf_surf_intersect_point<T, dim> current_point = current_intpoints.m_inter_points.back();
					Eigen::Vector2<T> next_param = intersect_param.template block<2, 1>(2, 0);
					Eigen::Vector2<T> test_param = boundary_point.m_uv.template block<2, 1>(2, 0);
					if (is_point_on_arc(current_point.m_uv, next_param, test_param, priori_enclosure, type, current_intpoints.m_is_positive_direction))
					{
						current_intpoint.m_priori_enclosure = priori_enclosure;
						// TODO: 处理priori_enclosure
						int points_count = m_boundary_chats[index].m_inter_points.size();
						for (int ipoint_index = points_count - 1; ipoint_index > 0; --ipoint_index)
						{
							m_boundary_chats[index].m_inter_points[ipoint_index].m_priori_enclosure = m_boundary_chats[index].m_inter_points[ipoint_index - 1].m_priori_enclosure;
						}

						current_intpoints.m_inter_points.insert(current_intpoints.m_inter_points.end(), m_boundary_chats[index].m_inter_points.rbegin(), m_boundary_chats[index].m_inter_points.rend());

						m_boundary_chats.erase(m_boundary_chats.begin() + index);
						is_stop = true;
						return ENUM_NURBS::NURBS_SUCCESS;
					}
				}
			}

			for (int index = 0; index < m_singular_chats.size(); ++index)
			{
				surf_surf_intersect_point<T, dim>& singular_point = m_singular_chats[index].m_inter_points.back();
				if (priori_enclosure.is_contain_point(singular_point.m_uv))
				{
					surf_surf_intersect_point<T, dim> current_point = current_intpoints.m_inter_points.back();
					Eigen::Vector2<T> next_param = intersect_param.template block<2, 1>(2, 0);
					Eigen::Vector2<T> test_param = singular_point.m_uv.template block<2, 1>(2, 0);
					if (is_point_on_arc(current_point.m_uv, next_param, test_param, priori_enclosure, type, current_intpoints.m_is_positive_direction))
					{
						current_intpoint.m_priori_enclosure = priori_enclosure;
						// TODO: 处理priori_enclosure
						int points_count = m_singular_chats[index].m_inter_points.size();
						for (int point_index = points_count - 1; point_index > 0; --point_index)
						{
							m_singular_chats[index].m_inter_points[point_index].m_priori_enclosure = m_singular_chats[index].m_inter_points[point_index - 1].m_priori_enclosure;
						}

						current_intpoints.m_inter_points.insert(current_intpoints.m_inter_points.end(), m_singular_chats[index].m_inter_points.rbegin(), m_singular_chats[index].m_inter_points.rend());
						m_singular_chats.erase(m_singular_chats.begin() + index);
						is_stop = true;
						return ENUM_NURBS::NURBS_SUCCESS;
					}
				}
			}

			surf_surf_intersect_point<T, dim> next_intpoint;
			next_intpoint.m_priori_enclosure = priori_enclosure;
			next_intpoint.m_uv = intersect_param;
			next_intpoint.m_previous_step = small_step;
			// m_left_ders(0, 0)->point_on_surface(intersect_param[0], intersect_param[1], next_intpoint.m_point);
			m_left_surf_compute->point_on_surface(m_left_ders(0, 0), intersect_param[0], intersect_param[1], next_intpoint.m_point);
			if (current_intpoints.start_point_type == POINT_TYPE::LOOP)
			{
				if (is_closed(current_intpoints, next_intpoint, type))
				{
					next_intpoint = current_intpoints.m_inter_points.front();
					current_intpoints.end_point_type = POINT_TYPE::LOOP;
					is_stop = true;
				}
			}

			// if (current_intpoints.start_point_type == POINT_TYPE::SINGULAR)
			// {
			// 	for (int index = 0; index < m_loop_chats.size(); ++index)
			// 	{
			// 		surf_surf_intersect_point<T, dim>& loop_point = m_loop_chats[index].m_inter_points.back();

			// 		if (priori_enclosure.is_contain_point(loop_point.m_uv))
			// 		{
			// 			surf_surf_intersect_point<T, dim> current_point = current_intpoints.m_inter_points.back();
			// 			Eigen::Vector2<T> next_param = intersect_param.template block<2, 1>(2, 0);
			// 			Eigen::Vector2<T> test_param = loop_point.m_uv.template block<2, 1>(2, 0);
			// 			if (is_point_on_arc(current_point.m_uv, next_param, test_param, priori_enclosure, type, current_intpoints.m_is_positive_direction))
			// 			{
			// 				current_intpoint.m_priori_enclosure = priori_enclosure;
			// 				// TODO: 处理priori_enclosure
			// 				int points_count = m_loop_chats[index].m_inter_points.size();
			// 				for (int ipoint_index = points_count - 1; ipoint_index > 0; --ipoint_index)
			// 				{
			// 					m_loop_chats[index].m_inter_points[ipoint_index].m_priori_enclosure = m_loop_chats[index].m_inter_points[ipoint_index - 1].m_priori_enclosure;
			// 				}
			// 				current_intpoints.m_inter_points.insert(current_intpoints.m_inter_points.end(), m_loop_chats[index].m_inter_points.rbegin(), m_loop_chats[index].m_inter_points.rend());

			// 				m_loop_chats.erase(m_loop_chats.begin() + index);
			// 				is_stop = true;
			// 				return ENUM_NURBS::NURBS_SUCCESS;
			// 			}
			// 		}
			// 	}
			// }
			current_intpoints.m_inter_points.back().m_priori_enclosure = priori_enclosure;
			current_intpoints.m_inter_points.push_back(next_intpoint);
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		ENUM_NURBS trace_curve_by_point(surfs_int_points_chat<T, dim>& current_intpoints)
		{
			int index = 0;
			bool is_closed;
			for (; index < MAXINTERSETORPOINTNUMBER; ++index)
			{
				if (trace_point(current_intpoints, is_closed) != ENUM_NURBS::NURBS_SUCCESS)
				{
					//TODO:析构内存
					return ENUM_NURBS::NURBS_ERROR;
				}
				if (is_closed == true)
				{
					break;
				}
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		ENUM_NURBS trace_int_chat()
		{
			ENUM_NURBS err_code;
			while (m_singular_chats.empty() == false)
			{
				surfs_int_points_chat<T, dim> current_intpoints = m_singular_chats.back();
				m_singular_chats.pop_back();
				err_code = trace_curve_by_point(current_intpoints);
				if (err_code != ENUM_NURBS::NURBS_SUCCESS)
				{
					return err_code;
				}
				m_result.m_int_chats.push_back(std::move(current_intpoints));
			}

			while (m_boundary_chats.empty() == false)
			{
				surfs_int_points_chat<T, dim> current_intpoints = m_boundary_chats.back();
				m_boundary_chats.pop_back();
				err_code = trace_curve_by_point(current_intpoints);
				if (err_code != ENUM_NURBS::NURBS_SUCCESS)
				{
					return err_code;
				}
				m_result.m_int_chats.push_back(std::move(current_intpoints));
			}

			while (m_loop_chats.empty() == false)
			{
				surfs_int_points_chat<T, dim> current_intpoints = m_loop_chats.back();
				m_loop_chats.pop_back();
				size_t index = 0;
				if (is_point_in_intcurve(current_intpoints.m_inter_points[0].m_uv, current_intpoints.m_is_transversal, index) == true)
				{
					continue;
				}
				else
				{
					err_code = trace_curve_by_point(current_intpoints);
					if (err_code != ENUM_NURBS::NURBS_SUCCESS)
					{
						return err_code;
					}
					m_result.m_int_chats.push_back(std::move(current_intpoints));
				}
			}

			return ENUM_NURBS::NURBS_SUCCESS;
		}

		ENUM_NURBS surf_surf_boundary_int(left_surface_type& lsurf, right_surface_type& rsurf,
			std::unordered_set<help<T, 4>, hash_help<T, 4>>& remove_mult)
		{
			right_curve curve1;
			right_curve curve2;
			right_curve curve3;
			right_curve curve4;
			std::array<T, 2> u_knots_end;
			std::array<T, 2> v_knots_end;
			rsurf.get_uv_knots_end(u_knots_end, v_knots_end);
			rsurf.template get_isoparameter_curve<ENUM_DIRECTION::U_DIRECTION>(u_knots_end[0], curve1);
			rsurf.template get_isoparameter_curve<ENUM_DIRECTION::U_DIRECTION>(u_knots_end[1], curve2);
			rsurf.template get_isoparameter_curve<ENUM_DIRECTION::V_DIRECTION>(v_knots_end[0], curve3);
			rsurf.template get_isoparameter_curve<ENUM_DIRECTION::V_DIRECTION>(v_knots_end[1], curve4);

			curve_surface_int<right_curve, left_surface_type> curve_surface_intersect;
			curve_surface_intersect.init(&curve1, &lsurf);
			curve_surface_intersect.run_intersect();

			for (auto& param : curve_surface_intersect.m_int_points)
			{
				remove_mult.insert(Eigen::Vector4<T>(param.m_int_param[1], param.m_int_param[2], u_knots_end[0], param.m_int_param[0]));
			}

			curve_surface_intersect.reset_curve(&curve2);
			curve_surface_intersect.run_intersect();
			for (auto& param : curve_surface_intersect.m_int_points)
			{
				remove_mult.insert(Eigen::Vector4<T>(param.m_int_param[1], param.m_int_param[2], u_knots_end[1], param.m_int_param[0]));
			}

			curve_surface_intersect.reset_curve(&curve3);
			curve_surface_intersect.run_intersect();
			for (auto& param : curve_surface_intersect.m_int_points)
			{
				remove_mult.insert(Eigen::Vector4<T>(param.m_int_param[1], param.m_int_param[2], param.m_int_param[0], v_knots_end[0]));
			}

			curve_surface_intersect.reset_curve(&curve4);
			curve_surface_intersect.run_intersect();
			for (auto& param : curve_surface_intersect.m_int_points)
			{
				remove_mult.insert(Eigen::Vector4<T>(param.m_int_param[1], param.m_int_param[2], param.m_int_param[0], u_knots_end[1]));
			}

			left_curve curve5;
			left_curve curve6;
			left_curve curve7;
			left_curve curve8;
			lsurf.get_uv_knots_end(u_knots_end, v_knots_end);
			lsurf.template get_isoparameter_curve<ENUM_DIRECTION::U_DIRECTION>(u_knots_end[0], curve5);
			lsurf.template get_isoparameter_curve<ENUM_DIRECTION::U_DIRECTION>(u_knots_end[1], curve6);
			lsurf.template get_isoparameter_curve<ENUM_DIRECTION::V_DIRECTION>(v_knots_end[0], curve7);
			lsurf.template get_isoparameter_curve<ENUM_DIRECTION::V_DIRECTION>(v_knots_end[1], curve8);
			curve_surface_int<left_curve, right_surface_type> curve_surface_intersect2;
			curve_surface_intersect2.init(&curve5, &rsurf);
			curve_surface_intersect2.run_intersect();

			for (auto& param : curve_surface_intersect2.m_int_points)
			{
				remove_mult.insert(Eigen::Vector4<T>(u_knots_end[0], param.m_int_param[0], param.m_int_param[1], param.m_int_param[2]));
			}

			curve_surface_intersect2.reset_curve(&curve6);
			curve_surface_intersect2.run_intersect();
			for (auto& param : curve_surface_intersect2.m_int_points)
			{
				remove_mult.insert(Eigen::Vector4<T>(u_knots_end[1], param.m_int_param[0], param.m_int_param[1], param.m_int_param[2]));
			}

			curve_surface_intersect2.reset_curve(&curve7);
			curve_surface_intersect2.run_intersect();
			for (auto& param : curve_surface_intersect2.m_int_points)
			{
				remove_mult.insert(Eigen::Vector4<T>(param.m_int_param[0], v_knots_end[0], param.m_int_param[1], param.m_int_param[2]));
			}

			curve_surface_intersect2.reset_curve(&curve8);
			curve_surface_intersect2.run_intersect();
			for (auto& param : curve_surface_intersect2.m_int_points)
			{
				remove_mult.insert(Eigen::Vector4<T>(param.m_int_param[0], u_knots_end[1], param.m_int_param[1], param.m_int_param[2]));
			}


			return ENUM_NURBS::NURBS_SUCCESS;
		}
		
		template<typename surface_type>
		ENUM_NURBS eval_ders_box(const Eigen::Vector2<T>& param, const surface_type& der10, const surface_type& der01, std::array<surface_type, 4>& der10_normal_boxes)
		{
			constexpr bool is_ratio = surface_type::is_ratio;
			std::array<surface_type, 4> der10_sub_surface, der01_sub_surface;
			der10->split_at_param(param[0], &der10_sub_surface);
			der01->split_at_param(param[1], &der01_sub_surface);
			return ENUM_NURBS::NURBS_SUCCESS;
		}


		template<typename surface_type>
		ENUM_NURBS get_2_box(const surface_type* surf, Box<T, dim>& box, Box<T, 2>& uv_box) const
		{
			surf->get_uv_box(uv_box);
			surf->get_box(box);
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		template<typename surface_type>
		ENUM_NURBS get_3_box(const surface_type* surf, Box<T, dim>& u_box, Box<T, dim>& v_box, Box<T, dim>& normal_box) const
		{
			if constexpr (surface_type::is_ratio == false)
			{
				surf->tangent_v_surface_box(v_box);
				surf->tangent_u_surface_box(u_box);
				normal_box = interval_algorithm::cross(u_box, v_box);
			}
			else
			{
				assert(false);
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		int surfs_patch_is_int(surfs_patch<left_surface_type>* surf1, surfs_patch<right_surface_type>* surf2)
		{
			// volatile bool f = surf1->m_box.intersect(surf2->m_box, intersect_box);
			if (surf1->m_box.intersect(surf2->m_box, intersect_box) == true)
			{
				if (surf1->m_box_is_valid == false)
				{
					get_3_box<left_surface_type>(surf1->m_surf, surf1->m_u_tangent_box, surf1->m_v_tangent_box, surf1->m_normal_box);
					surf1->m_box_is_valid = true;
				}
				if (surf2->m_box_is_valid == false)
				{
					get_3_box<right_surface_type>(surf2->m_surf, surf2->m_u_tangent_box, surf2->m_v_tangent_box, surf2->m_normal_box);
					surf2->m_box_is_valid = true;
				}
				Box<T, 4> param_ders;

				Box<T, 3> space_ders;
				bool may_arrived_singular = false;
				Box<T, 4> param_ders_box;
				ENUM_NURBS code = eval_preiamge_and_space_ders(surf1->m_u_tangent_box, surf1->m_v_tangent_box, surf2->m_u_tangent_box, surf2->m_v_tangent_box,
					surf1->m_normal_box, surf2->m_normal_box, space_ders, param_ders, may_arrived_singular);
				if (code != ENUM_NURBS::NURBS_SUCCESS)
				{
					return 2;
				}
				Box<T, 2> uv_tangent_box(param_ders.Min.template block<2, 1>(0, 0), param_ders.Max.template block<2, 1>(0, 0));

				if (uv_tangent_box.is_contain_point(Eigen::Vector2<T>(0, 0)))
				{
					return 1;
				}
				// cotian(-1. 0)
				// if (uv_tangent_box.Min[1] <= PRECISION<T>::value && uv_tangent_box.Max[1] >= -PRECISION<T>::value && uv_tangent_box.Min[0] <= PRECISION<T>::value)
				// {
				// 	return 1;
				// }
				// cotian(0. 1)
				if (uv_tangent_box.Min[0] <= PRECISION<T>::value && uv_tangent_box.Max[0] >= -PRECISION<T>::value && uv_tangent_box.Min[1] >= -PRECISION<T>::value)
				{
					return 1;
				}
				// std::string stringx = "surface" + std::to_string(surf_index) + ".obj";
				// const char* path = stringx.c_str();
				// save_obj2(*surf1->m_surf, path);
				// surf_index += 1;
				// std::string stringx2 = "surface" + std::to_string(surf_index) + ".obj";
				// const char* path2 = stringx2.c_str();
				// save_obj2(*surf2->m_surf, path2);
				// surf_index += 1;
				return 0;
			}
			return 0;
		}


		ENUM_NURBS get_initial_box(std::vector<std::pair<Box<T, 2>, Box<T, 2>>>& sub_int_boxes, std::vector<char>& is_tangent)
		{
			sub_int_boxes.clear();
			is_tangent.clear();
			
			std::vector<left_surface_type> left_surface_int;
			std::vector<right_surface_type> right_surface_int;

			Box<T, dim> left_sub_boxes, right_sub_boxes;
			Box<T, dim> left_sub_normal_boxes, right_sub_normal_boxes;
			Box<T, dim> u_left_sub_tangent_boxes, u_right_sub_tangent_boxes;
			Box<T, dim> v_left_sub_tangent_boxes, v_right_sub_tangent_boxes;
			std::vector<left_surface_type*> left_sub_surfaces;
			std::vector<right_surface_type*> right_sub_surfaces;
			Box<T, 2> uv_box, st_box;

			left_surface_type* left_surface_temp = new left_surface_type(*m_left_ders(0, 0));
			right_surface_type* right_surface_temp = new right_surface_type(*m_right_ders(0, 0));

			surfs_patch<left_surface_type>* left_surf_patch_root = new surfs_patch<left_surface_type>();
			surfs_patch<right_surface_type>* right_surf_patch_root = new surfs_patch<right_surface_type>();
			left_surf_patch_root->m_surf = left_surface_temp;
			right_surf_patch_root->m_surf = right_surface_temp;
			get_2_box<left_surface_type>(left_surface_temp, left_surf_patch_root->m_box, left_surf_patch_root->m_uv_box);
			get_2_box<right_surface_type>(right_surface_temp, right_surf_patch_root->m_box, right_surf_patch_root->m_uv_box);
			
			int_surfs_pair<left_surface_type, right_surface_type>* int_pair = new int_surfs_pair<left_surface_type, right_surface_type>();
			int_pair->m_left_surfs_patch = left_surf_patch_root;
			int_pair->m_right_surfs_patch = right_surf_patch_root;

			Eigen::Vector2<T> uv_param;
			Eigen::Vector2<T> st_param;
			Box<T, dim> intersect_box;

			std::vector<surfs_patch<left_surface_type>*> current_int_patches;
			if (surfs_patch_is_int(left_surf_patch_root, right_surf_patch_root) != 0)
			{
				left_surf_patch_root->m_int_patches.push_back(right_surf_patch_root);
				current_int_patches.push_back(left_surf_patch_root);
			}
			std::vector<surfs_patch<left_surface_type>*> next_int_patches;
			int deep = 0;
			int max_deep = 8;
			while (current_int_patches.size() != 0)
			{
				++deep;
				for (surfs_patch<left_surface_type>* left_int_surf : current_int_patches)
				{
					uv_param = left_int_surf->m_uv_box.get_middle_point();
					// left_int_surf->m_surf->split_at_param(uv_param, left_sub_surfaces);
					m_left_surf_compute->split_at_param(left_int_surf->m_surf, uv_param, left_sub_surfaces);
					left_int_surf->m_children = std::vector<surfs_patch<left_surface_type>*>{ new surfs_patch<left_surface_type>{}, new surfs_patch<left_surface_type>{}, new surfs_patch<left_surface_type>{}, new surfs_patch<left_surface_type>{} };
					for (int index = 0; index < 4; ++index)
					{
						left_int_surf->m_children[index]->m_surf = left_sub_surfaces[index];
						get_2_box<left_surface_type>(left_int_surf->m_children[index]->m_surf, left_int_surf->m_children[index]->m_box, left_int_surf->m_children[index]->m_uv_box);
					}

					std::vector<surfs_patch<right_surface_type>*> right_int_surfs = left_int_surf->m_int_patches;
					for (surfs_patch<right_surface_type>* right_int_surf : right_int_surfs)
					{
						if (right_int_surf->m_children.size() == 0)
						{
							st_param = right_int_surf->m_uv_box.get_middle_point();
							// right_int_surf->m_surf->split_at_param(st_param, right_sub_surfaces);
							m_right_surf_compute->split_at_param(right_int_surf->m_surf, st_param, right_sub_surfaces);
							right_int_surf->m_children = std::vector<surfs_patch<right_surface_type>*>{ new surfs_patch<right_surface_type>(), new surfs_patch<right_surface_type>(), new surfs_patch<right_surface_type>(), new surfs_patch<right_surface_type>() };
							for (int index = 0; index < 4; ++index)
							{
								right_int_surf->m_children[index]->m_surf = right_sub_surfaces[index];
								get_2_box<right_surface_type>(right_int_surf->m_children[index]->m_surf, right_int_surf->m_children[index]->m_box, right_int_surf->m_children[index]->m_uv_box);
							}
						}
						for (surfs_patch<left_surface_type>* left_int_surf_children : left_int_surf->m_children)
						{
							for (surfs_patch<right_surface_type>* right_int_surf_children : right_int_surf->m_children)
							{
								int int_type = surfs_patch_is_int(left_int_surf_children, right_int_surf_children);
								if (deep == max_deep)
								{
									if (int_type == 1)
									{
										sub_int_boxes.push_back(std::make_pair(left_int_surf_children->m_uv_box, right_int_surf_children->m_uv_box));
										is_tangent.push_back('0');
									}
									else if (int_type == 2)
									{
										sub_int_boxes.push_back(std::make_pair(left_int_surf_children->m_uv_box, right_int_surf_children->m_uv_box));
										is_tangent.push_back('1');
									}
								}
								else
								{
									if (int_type != 0)
									{
										left_int_surf_children->m_int_patches.push_back(right_int_surf_children);
									}
								}
							}
						}
					}
					if (deep < max_deep)
					{
						for (surfs_patch<left_surface_type>* left_int_surf_children : left_int_surf->m_children)
						{
							if (left_int_surf_children->m_int_patches.size() > 0)
							{
								next_int_patches.push_back(left_int_surf_children);
							}
						}
					}
				}
				std::swap(current_int_patches, next_int_patches);
				next_int_patches.clear();
			}
			
			delete int_pair;
			return ENUM_NURBS::NURBS_SUCCESS;

		};
		

		ENUM_NURBS get_initial_box2(std::vector<std::pair<Box<T, 2>, Box<T, 2>>>& sub_int_boxes, std::vector<char>& is_tangent)
		{
			sub_int_boxes.clear();
			is_tangent.clear();
			
			// std::vector<std::pair<left_surface_type, right_surface_type>> int_surface_pairs;
		
			std::vector<left_surface_type> left_surface_int;
			std::vector<right_surface_type> right_surface_int;
			// std::vector<left_surface_type_sptr> left_surface_int;
			// std::vector<right_surface_type_sptr> right_surface_int;

			Box<T, dim> left_sub_boxes, right_sub_boxes;
			Box<T, dim> left_sub_normal_boxes, right_sub_normal_boxes;
			Box<T, dim> u_left_sub_tangent_boxes, u_right_sub_tangent_boxes;
			Box<T, dim> v_left_sub_tangent_boxes, v_right_sub_tangent_boxes;
			std::array<left_surface_type, 4> left_sub_surfaces;
			std::array<right_surface_type, 4> right_sub_surfaces;
			// std::array<left_surface_type_sptr, 4> left_sub_surfaces{ nullptr, nullptr, nullptr, nullptr };
			// std::array<right_surface_type_sptr, 4> right_sub_surfaces{ nullptr, nullptr, nullptr, nullptr };
			Box<T, 2> uv_box, st_box;
		    m_left_ders(0, 0)->get_uv_box(uv_box);
			m_right_ders(0, 0)->get_uv_box(st_box);
			Eigen::Vector2<T> uv_param = uv_box.get_middle_point();
			Eigen::Vector2<T> st_param = st_box.get_middle_point();

			left_surface_type left_surface_temp(*m_left_ders(0, 0));
			right_surface_type right_surface_temp(*m_right_ders(0, 0));

			// m_left_ders(0, 0)->split_at_param(uv_param, left_sub_surfaces);
			// m_right_ders(0, 0)->split_at_param(st_param, right_sub_surfaces);
			
			// std::map<left_surface_type_sptr, std::array<left_surface_type_sptr, 4>> left_pairs;
			// std::map<right_surface_type_sptr, std::array<right_surface_type_sptr, 4>> right_pairs;

			// left_surface_temp.split_at_param(uv_param, left_sub_surfaces2);
			left_surface_temp.split_at_param(uv_param, left_sub_surfaces);
			// left_pairs[left_surface_temp] = left_sub_surfaces;
			
			right_surface_temp.split_at_param(st_param, right_sub_surfaces);
			// right_pairs[right_surface_temp] = right_sub_surfaces;
			Box<T, dim> intersect_box;
			for (int i = 0; i < 4; ++i)
			{
				// left_sub_surfaces[i]->get_box(left_sub_boxes);
				left_sub_surfaces[i].get_box(left_sub_boxes);
				for (int j = 0; j < 4; ++j)
				{
					// right_sub_surfaces[j]->get_box(right_sub_boxes);
					right_sub_surfaces[j].get_box(right_sub_boxes);
					if (left_sub_boxes.intersect(right_sub_boxes, intersect_box))
					{
						left_surface_int.push_back(left_sub_surfaces[i]);
						right_surface_int.push_back(right_sub_surfaces[j]);
						// int_surface_pairs.push_back(std::make_pair(left_sub_surfaces[i], right_sub_surfaces[j]));
					}
				}
			}
			// while (int_surface_pairs.empty() == false)
			while (left_surface_int.empty() == false)
			{
				// std::pair<left_surface_type, right_surface_type> surface_pair = int_surface_pairs.back();
				// int_surface_pairs.pop_back();
				
				left_surface_type& left_surf = left_surface_int.back();
				right_surface_type& right_surf = right_surface_int.back();
				// left_surface_type_sptr left_surf = left_surface_int.back();
				// right_surface_type_sptr right_surf = right_surface_int.back();

				// left_surf->get_uv_box(uv_box);
				// right_surf->get_uv_box(st_box);
				left_surf.get_uv_box(uv_box);
				right_surf.get_uv_box(st_box);

				// surface_pair.first.get_uv_box(uv_box);
				// surface_pair.second.get_uv_box(st_box);
				Eigen::Vector2<T> uv_middle = uv_box.get_middle_point();
				Eigen::Vector2<T> st_middle = st_box.get_middle_point();
				// surface_pair.first.split_at_param(uv_middle, left_sub_surfaces);
				// surface_pair.second.split_at_param(st_middle, right_sub_surfaces);
				// left_surf->split_at_param(uv_middle, left_sub_surfaces);
				// right_surf->split_at_param(st_middle, right_sub_surfaces);
				left_surf.split_at_param(uv_middle, left_sub_surfaces);
				right_surf.split_at_param(st_middle, right_sub_surfaces);
				// if (left_surf->is_valid == false)
				// {
				// 	left_sub_surfaces = left_pairs[left_surf];
				// }
				// else
				// {
				// 	left_surf->split_at_param(uv_middle, left_sub_surfaces);
				// 	// left_pairs[left_surf] = left_sub_surfaces;
				// }
				// 
				// if (right_surf->is_valid == false)
				// {
				// 	right_sub_surfaces = right_pairs[right_surf];
				// }
				// else
				// {
				// 	right_surf->split_at_param(st_middle, right_sub_surfaces);
				// 	// right_pairs[right_surf] = right_sub_surfaces;
				// }

				left_surface_int.pop_back();
				right_surface_int.pop_back();

				T temp = 1.0 / 2.0 + 1e-4;
				for (int i = 0; i < 4; ++i)
				{
					left_sub_surfaces[i].get_box(left_sub_boxes);
					left_sub_surfaces[i].get_uv_box(uv_box);
					// left_sub_surfaces[i]->get_box(left_sub_boxes);
					// left_sub_surfaces[i]->get_uv_box(uv_box);
					for (int j = 0; j < 4; ++j)
					{
						// right_sub_surfaces[j]->get_uv_box(st_box);
						// right_sub_surfaces[j]->get_box(right_sub_boxes);
						right_sub_surfaces[j].get_uv_box(st_box);
						right_sub_surfaces[j].get_box(right_sub_boxes);
						if (left_sub_boxes.intersect(right_sub_boxes, intersect_box) == true)
						{
							Box<T, 4> param_ders;

							Box<T, 3> space_ders;
							bool may_arrived_singular = false;
							Box<T, 4> param_ders_box;

							// eval_normal_and_ders_box<left_surface_type>(*left_sub_surfaces[i], u_left_sub_tangent_boxes, u_right_sub_tangent_boxes, left_sub_normal_boxes);
							// eval_normal_and_ders_box<right_surface_type>(*right_sub_surfaces[j], v_left_sub_tangent_boxes, v_right_sub_tangent_boxes, right_sub_normal_boxes);
							eval_normal_and_ders_box<left_surface_type>(left_sub_surfaces[i], u_left_sub_tangent_boxes, v_left_sub_tangent_boxes, left_sub_normal_boxes);
							eval_normal_and_ders_box<right_surface_type>(right_sub_surfaces[j], u_right_sub_tangent_boxes, v_right_sub_tangent_boxes, right_sub_normal_boxes);
							ENUM_NURBS code = eval_preiamge_and_space_ders(u_left_sub_tangent_boxes, v_left_sub_tangent_boxes, u_right_sub_tangent_boxes, v_right_sub_tangent_boxes, left_sub_normal_boxes, right_sub_normal_boxes, space_ders, param_ders, may_arrived_singular);
							if (code != ENUM_NURBS::NURBS_SUCCESS)
							{
								// 目前对相切的时候不做处理; TODO: 找一下优化方法？
								if (std::abs(uv_box.Min[0] - uv_box.Max[0]) < temp)
								{
									bool add_int_box = true;
									Box<T, 2> intersect_box1, intersect_box2;
									for (Box<T, 4> &singular_box : m_singular_boxes)
									{
										Box<T, 2> uv_singular_box(singular_box.Min.template block<2, 1>(0, 0), singular_box.Max.template block<2, 1>(0, 0));
										Box<T, 2> st_singular_box(singular_box.Min.template block<2, 1>(2, 0), singular_box.Max.template block<2, 1>(2, 0));
										if (uv_box.intersect(uv_singular_box, intersect_box1) || st_box.intersect(st_singular_box, intersect_box2))
										{
											add_int_box = false;
											break;
										}
									}
									if (add_int_box == true)
									{
										sub_int_boxes.push_back(std::make_pair(uv_box, st_box));
										is_tangent.push_back('1');
									}

								}
								else
								{
									left_surface_int.push_back(left_sub_surfaces[i]);
									right_surface_int.push_back(right_sub_surfaces[j]);
									// int_surface_pairs.push_back(std::make_pair(left_sub_surfaces[i], right_sub_surfaces[j]));
									
								}
								continue;
							}

							Box<T, 2> uv_tangent_box(param_ders.Min.template block<2, 1>(0, 0), param_ders.Max.template block<2, 1>(0, 0));
							// cone<T, 2> uv_tangent_cone = point_box(uv_tangent_box, Eigen::Vector2<T>(0, 0));

							// cotian(-1. 0)
							if (uv_tangent_box.Min[1] <= PRECISION<T>::value && uv_tangent_box.Max[1] >= -PRECISION<T>::value && uv_tangent_box.Min[0] <= PRECISION<T>::value)
							{
								if (std::abs(uv_box.Min[0] - uv_box.Max[0]) < temp)
								{
									bool add_int_box = true;

									Box<T, 2> intersect_box1, intersect_box2;
									for (Box<T, 4> &singular_box : m_singular_boxes)
									{
										Box<T, 2> uv_singular_box(singular_box.Min.template block<2, 1>(0, 0), singular_box.Max.template block<2, 1>(0, 0));
										Box<T, 2> st_singular_box(singular_box.Min.template block<2, 1>(2, 0), singular_box.Max.template block<2, 1>(2, 0));
										if (uv_box.intersect(uv_singular_box, intersect_box1) || st_box.intersect(st_singular_box, intersect_box2))
										{
											add_int_box = false;
											break;
										}
									}

									if (add_int_box)
									{
										sub_int_boxes.push_back(std::make_pair(uv_box, st_box));
										is_tangent.push_back('0');
									}
								}
								else
								{
									left_surface_int.push_back(left_sub_surfaces[i]);
									right_surface_int.push_back(right_sub_surfaces[j]);
									// int_surface_pairs.push_back(std::make_pair(left_sub_surfaces[i], right_sub_surfaces[j]));
								}
							}
						}
					}
				}
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		};
		

		ENUM_NURBS find_inner_intersect_points(const std::vector<std::pair<Box<T, 2>, Box<T, 2>>>& sub_int_boxes, const std::vector<char>& is_tangent)
		{
			size_t intersect_count = sub_int_boxes.size();
			for (size_t index = 0; index < intersect_count; ++index)
			{
				const std::pair<Box<T, 2>, Box<T, 2>>& int_box = sub_int_boxes[index];
				Eigen::Vector4<T> initial_param, intersect_param;
				initial_param.template block<2, 1>(0, 0) = int_box.first.get_middle_point();
				initial_param.template block<2, 1>(2, 0) = int_box.second.get_middle_point();
				int type = 0;
				if (intersect_point_iteration(m_product_box, initial_param, intersect_param) == ENUM_NURBS::NURBS_SUCCESS)
				{
					if (false == int_box.first.is_contain_point(intersect_param.template block<2, 1>(0, 0))
						|| false == int_box.second.is_contain_point(intersect_param.template block<2, 1>(2, 0)))
					{
						continue;
					}
					if (is_tangent[index] == '1')
					{
						Eigen::Vector<T, 4> singular_intersect_param;
						if (ENUM_NURBS::NURBS_SUCCESS == intersect_singular_point_iteration(m_product_box, intersect_param, singular_intersect_param))
						{
							// if (int_box.first.is_contain_point(singular_intersect_param.template block<2, 1>(0, 0), PRECISION<T>::value) == false
							//     || int_box.second.is_contain_point(singular_intersect_param.template block<2, 1>(2, 0), PRECISION<T>::value) == false)
							// {
							//     continue;
							// }

							std::vector<Eigen::Matrix<T, 4, 1>> param_ders_temp;
							std::vector<Eigen::Matrix<T, dim, 1>> space_ders_temp;
							std::array<T, 4> ks;
							std::array<T, 2> left_EG, right_EG;
							ENUM_NURBS error_code = eval_preiamge_and_space_ders<1>(singular_intersect_param, param_ders_temp, space_ders_temp, ks, left_EG, right_EG, 3);
							if (error_code == ENUM_NURBS::NURBS_ISOLATED_TANGENTIAL_POINT)
							{
								surf_surf_intersect_point<T, dim> isolate_point;
								isolate_point.m_uv = singular_intersect_param;
								// m_right_ders(0, 0)->point_on_surface(singular_intersect_param[2], singular_intersect_param[3], isolate_point.m_point);
								m_right_surf_compute->point_on_surface(m_right_ders(0, 0), singular_intersect_param[2], singular_intersect_param[3], isolate_point.m_point);
								m_result.m_isolate_points.push_back(isolate_point);
								continue;
							}
							else if (param_ders_temp.size() == 1)
							{
								surfs_int_points_chat<T, dim> loop_points_chat;
								loop_points_chat.start_point_type = loop_points_chat.end_point_type = POINT_TYPE::LOOP;
								loop_points_chat.m_is_transversal = false;
								loop_points_chat.m_is_positive_direction = true;
								std::vector<surf_surf_intersect_point<T, dim>> loop_point(1);
								loop_point[0].m_uv = singular_intersect_param;
								// m_left_ders(0, 0)->point_on_surface(singular_intersect_param[0], singular_intersect_param[1], loop_point[0].m_point);
								m_left_surf_compute->point_on_surface(m_left_ders(0, 0), singular_intersect_param[0], singular_intersect_param[1], loop_point[0].m_point);
								loop_points_chat.m_inter_points = loop_point;
								m_loop_chats.push_back(std::move(loop_points_chat));
							}
							else
							{
								eval_singular_points(singular_intersect_param, param_ders_temp, space_ders_temp, m_singular_chats);
							}
						}
					}
					else
					{
						surfs_int_points_chat<T, dim> loop_points_chat;
						loop_points_chat.start_point_type = loop_points_chat.end_point_type = POINT_TYPE::LOOP;
						loop_points_chat.m_is_transversal = true;
						loop_points_chat.m_is_positive_direction = true;
						std::vector<surf_surf_intersect_point<T, dim>> loop_point(1);
						loop_point[0].m_uv = intersect_param;
						// m_left_ders(0, 0)->point_on_surface(intersect_param[0], intersect_param[1], loop_point[0].m_point);
						m_left_surf_compute->point_on_surface(m_left_ders(0, 0), intersect_param[0], intersect_param[1], loop_point[0].m_point);
						loop_points_chat.m_inter_points = loop_point;
						m_loop_chats.push_back(std::move(loop_points_chat));

					}
				}
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		// for non-singular points
		ENUM_NURBS make_int_chat_by_param(const std::unordered_set<help<T, 4>, hash_help<T, 4>>& params, std::vector<surfs_int_points_chat<T, dim>>& int_chats)
		{
			for (auto& param : params)
			{
				Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> left_ders;
				// m_left_ders(0, 0)->template derivative_on_surface<1>(param.m_point[0], param.m_point[1], left_ders);
				m_left_surf_compute->template derivative_on_surface<1>(param.m_point[0], param.m_point[1], left_ders);
				Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> right_ders;
				// m_right_ders(0, 0)->template derivative_on_surface<1>(param.m_point[2], param.m_point[3], right_ders);
				m_right_surf_compute->template derivative_on_surface<1>(param.m_point[2], param.m_point[3], right_ders);
				std::vector<surf_surf_intersect_point<T, dim>> current_intpoint(1);
				current_intpoint[0].m_uv = param.m_point;
				current_intpoint[0].m_point = left_ders(0, 0);
				Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
				left_normal.normalize();
				Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
				right_normal.normalize();

				bool transversal = true;
				double cosTheta = left_normal.dot(right_normal);
				double theta = std::acos(cosTheta);
				if (std::abs(cosTheta) > 1.0 - KNOTS_EPS<T>::value)
				{
					transversal = false;
				}
				surfs_int_points_chat<T, dim> boundary_points;
				boundary_points.m_inter_points = current_intpoint;
				boundary_points.m_is_transversal = transversal;
				boundary_points.m_is_positive_direction = true;
				boundary_points.start_point_type = POINT_TYPE::BOUNDARY;
				bool is_closed = false;
				if (trace_point(boundary_points, is_closed) == ENUM_NURBS::NURBS_SUCCESS)
				{
					int_chats.push_back(boundary_points);
				}
				assert(is_closed == false);
				boundary_points.m_is_positive_direction = false;
				boundary_points.m_inter_points = current_intpoint;
				if (trace_point(boundary_points, is_closed) == ENUM_NURBS::NURBS_SUCCESS)
				{
					int_chats.push_back(boundary_points);
				}
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		//1. 目前如果奇异点所连接的一条交线是相切交线，可能会丢失奇异点。导致丢失交线。应该可以使用细分的找到奇异点
		//2. 求相切的交线太耗时
		ENUM_NURBS surafces_intersection2()
		{

#ifdef DSSI
			using namespace std::chrono;
			auto s1 = steady_clock::now();
#endif // DSSI

			std::unordered_set<help<T, 4>, hash_help<T, 4>> remove_mult;
			surf_surf_boundary_int(*m_left_ders(0, 0), *m_right_ders(0, 0), remove_mult);
#ifdef DSSI
			auto e1 = steady_clock::now();
			auto time1 = duration_cast<microseconds>(e1 - s1);
			std::cout << "surf_surf_boundary_int: " <<  time1.count() << "um" << std::endl;
			auto s2 = steady_clock::now();
#endif // DSSI
			make_int_chat_by_param(remove_mult, m_boundary_chats);
#ifdef DSSI
			auto e2 = steady_clock::now();
			auto time2 = duration_cast<microseconds>(e2 - s2);
			std::cout << "make_int_chat_by_param: " <<  time2.count() << "um" << std::endl;
			auto s3 = steady_clock::now();
#endif // DSSI
			std::vector<std::pair<Box<T, 2>, Box<T, 2>>> sub_int_boxes;
			std::vector<char> is_tangent;
			get_initial_box(sub_int_boxes, is_tangent);

#ifdef DSSI
			auto e3 = steady_clock::now();
			auto time3 = duration_cast<microseconds>(e3 - s3);
			std::cout << "get_initial_box: " << time3.count() << "um" << std::endl;
			for (auto& int_pair : sub_int_boxes)
			{
				Box<T, 2> b = int_pair.first;
				Eigen::Vector2<T> uv_param = b.get_middle_point();
				Eigen::Vector3<T> p;
				m_left_ders(0, 0)->point_on_surface(uv_param[0], uv_param[1], p);
				loop_points.push_back(p);
			}

			auto s4 = steady_clock::now();
#endif // DSSI
			find_inner_intersect_points(sub_int_boxes, is_tangent);
#ifdef DSSI
			auto e4 = steady_clock::now();
			auto time4 = duration_cast<microseconds>(e4 - s4);
			std::cout << "find_inner_intersect_points: " << time4.count() << "um" << std::endl;
			auto s5 = steady_clock::now();
#endif // DSSI
			trace_int_chat();
#ifdef DSSI
			auto e5 = steady_clock::now();
			std::cout << "trace_int_chat: " << duration_cast<microseconds>(e5 - s5).count() << std::endl;
#endif // DSSI
			return ENUM_NURBS::NURBS_SUCCESS;
		}
	};
}







