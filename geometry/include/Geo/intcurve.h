#pragma once


#include "declare.h"
#include "nurbs_surface.h"
#include "interval_alogrithm.h"
#include <memory>
#include "bezier_curve_int.h"
#include <ctime>
#include <chrono>
#include <vector>
#include "vodes.h"
#include "surface.h"
#include <array>
#include <algorithm>
#include <type_traits>
#include <unsupported/Eigen/Polynomials>

namespace tnurbs
{

	template<typename Scalar>
	struct vec3_utils
	{
		using vec_type = Eigen::Vector3<Scalar>;
		// ͨ��һ����λ����(���뱣֤)���ҵ������һ���໥�����ĵ�λ������ʹ�ù���һ�鵥λ������
		// basis[0] = x_dir, basis[1] = y_dir
		static ENUM_NURBS create_lcs(const vec_type& normal_vec, std::array<vec_type, 2>& basis)
		{
			if (std::abs(normal_vec[0]) > 0.5)
			{
				basis[0] = vec_type(-normal_vec[1], normal_vec[0], 0.0);
			}
			else
			{
				basis[0] = vec_type(0.0, -normal_vec[2], normal_vec[1]);
			}
			basis[0].normalize();
			basis[1] = normal_vec.cross(basis[0]);
			basis[1].normalize();
			return ENUM_NURBS::NURBS_SUCCESS;
		}
	};



	template<typename Scalar>
	struct compare
	{
		static bool is_equal(const Scalar a, const Scalar b, const Scalar tol)
		{
			Scalar c = a - b;
			if (c > tol || c < -tol)
			{
				return false;
			}
			return true;
		}
		// a < b + tol ?
		static bool is_less(const Scalar a, const Scalar b, const Scalar tol)
		{
			if (a < b + tol)
			{
				return true;
			}
			return false;
		}
		
		//a > b + tol
		static bool is_greater(const Scalar a, const Scalar b, const Scalar tol)
		{
			if (a > b + tol)
			{
				return true;
			}
			return false;
		}
	};
	template<typename Scalar>
	struct solve_plolynomial_2_degree
	{
		// b11 * x^2 + b12 * x + b13
		static std::vector<Scalar> solve(const Scalar b11, const Scalar b12, const Scalar b22)
		{
			Scalar delta = b11 * b22 - b12 * b12;
			std::vector<Scalar> real_roots;
			auto check_result = [&](const Scalar real_root) -> bool
				{
					Scalar result = real_root * real_root * b11 * b11 + 2.0 * b12 * real_root + b22;
					if (compare<Scalar>::is_equal(result, 0.0, b11 * b11 * PRECISION<Scalar>::value * 1e-3))
					{
						return true;
					}
					return false;
				};

			if (compare<Scalar>::is_equal(delta, 0, b11 * b11 * PRECISION<Scalar>::value * 1e-3))
			{
				Scalar real_root = -b12 / (2.0 * b11);

				if (std::isnan(real_root))
				{
					assert(false);
					return ENUM_NURBS::NURBS_ERROR;
				}
				if (check_result(real_root) == false)
				{
					return ENUM_NURBS::NURBS_ERROR;
				}
				real_roots.push_back(real_root);
			}
			else
			{
				if (compare<Scalar>::is_equal(delta, 0, b11 * b11 * PRECISION<Scalar>::value))
				{
					// eigen ����ʽ���
					Eigen::Vector3<Scalar> coeff(b11, b12, b22);
					Eigen::PolynomialSolver<Scalar, 2> solver;
					solver.compute();
					const auto& roots = solver.roots();
					for (const auto& root : roots)
					{
						// check real result
						Scalar real_root = root.real();
						if (check_result(real_root) == true)
						{
							real_roots.push_back(real_root);
						}
					}
				}
			}

			return real_roots;
		}
	};

	template<typename T, int dim>
	struct int_curve_evaluate
	{
		using surf_type = surface<T, dim>;
		// tol must be positive
		static ENUM_NURBS evaluate_tangent_normal(T u, T v, T s, T t, const surf_type& left_surf, const surf_type& right_surf, Eigen::Vector<T, dim>& tangent)
		{
			// Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> left_ders;
			// ENUM_NURBS err_code = left_surf.derivative_on_surface<1>(u, v, left_ders);
			Eigen::MatrixX<Eigen::Vector<T, dim>> left_ders;
			ENUM_NURBS err_code = left_surf.derivative_on_surface(1, u, v, left_ders);
			if (err_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				return err_code;
			}
			// Eigen::Matrix<Eigen::Vector<T, dim>, 2, 2> right_ders;
			// err_code = right_surf.derivative_on_surface<1>(s, t, right_ders);
			Eigen::MatrixX<Eigen::Vector<T, dim>> right_ders;
			err_code = right_surf.derivative_on_surface(1, s, t, right_ders);
			if (err_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				return err_code;
			}

			// normal direction
			Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
			if (left_normal.squaredNorm() < PRECISION<T>::value * 1e-6)
			{
				return ENUM_NURBS::NURBS_DEGENERATE;
			}
			left_normal.normalize();

			Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
			if (right_normal.squaredNorm() < PRECISION<T>::value * 1e-6)
			{
				return ENUM_NURBS::NURBS_DEGENERATE;
			}
			right_normal.normalize();

			tangent = left_normal.cross(right_normal);
			if (tangent.squaredNorm() < PRECISION<T>::value * 1e-6)
			{
				return ENUM_NURBS::NURBS_DEGENERATE;
			}
			tangent.normalize();
			return ENUM_NURBS::NURBS_DEGENERATE;
		}

		static ENUM_NURBS evaluate_tangent_parellel(T u, T v, T s, T t, const surf_type& left_surf, const surf_type& right_surf, std::vector<Eigen::Vector<T, dim>>& tangents)
		{

			Eigen::MatrixX<Eigen::Vector<T, dim>> left_ders;
			ENUM_NURBS err_code = left_surf.derivative_on_surface(2, u, v, left_ders);
			// Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> left_ders;
			// ENUM_NURBS err_code = left_surf.derivative_on_surface<2>(u, v, left_ders);
			if (err_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				return err_code;
			}
			// Eigen::Matrix<Eigen::Vector<T, dim>, 3, 3> right_ders;
			// err_code = right_surf.derivative_on_surface<3>(u, v, right_ders);
			Eigen::MatrixX<Eigen::Vector<T, dim>> right_ders;
			err_code = right_surf.derivative_on_surface(2, u, v, right_ders);
			if (err_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				return err_code;
			}

			// normal direction
			Eigen::Vector<T, dim> left_normal = left_ders(1, 0).cross(left_ders(0, 1));
			if (left_normal.squaredNorm() < PRECISION<T>::value * 1e-6)
			{
				return ENUM_NURBS::NURBS_DEGENERATE;
			}
			left_normal.normalize();

			Eigen::Vector<T, dim> right_normal = right_ders(1, 0).cross(right_ders(0, 1));
			if (right_normal.squaredNorm() < PRECISION<T>::value * 1e-6)
			{
				return ENUM_NURBS::NURBS_DEGENERATE;
			}
			right_normal.normalize();

			T vec_dot = left_normal.dot(right_normal);
			if (std::abs(vec_dot - 1.0) > PRECISION<T>::value * 1e-6)
			{
				return ENUM_NURBS::NURBS_NO_TANGENT;
			}
			// bool is_same_dir = (vec_dot > 0.0);

			T Eb = right_ders(1, 0).dot(right_ders(1, 0));
			T Fb = right_ders(1, 0).dot(right_ders(0, 1));
			T Gb = right_ders(0, 1).dot(right_ders(0, 1));

			T denominator = std::sqrt(Eb * Gb - Fb * Fb);
			if (std::isnan(denominator))
			{
				assert(false);
				return ENUM_NURBS::NURBS_ERROR;
			}
			if (is_equal(denominator, 0.0, PRECISION<T>::value))
			{
				return ENUM_NURBS::NURBS_ERROR;
			}

			T denominator_r = 1.0 / denominator;
			T a11 = (left_ders(1, 0).cross(right_ders(0, 1))).dot(left_normal) * denominator_r;
			T a12 = (left_ders(0, 1).cross(right_ders(0, 1))).dot(left_normal) * denominator_r;
			T a21 = (right_ders(1, 0).cross(left_ders(1, 0))).dot(left_normal) * denominator_r;
			T a22 = (right_ders(1, 0).cross(left_ders(0, 1))).dot(left_normal) * denominator_r;

			T La = left_ders(2, 0).dot(left_normal);
			T Ma = left_ders(1, 1).dot(left_normal);
			T Na = left_ders(0, 2).dot(left_normal);

			T Lb = right_ders(2, 0).dot(left_normal);
			T Mb = right_ders(1, 1).dot(left_normal);
			T Nb = right_ders(0, 2).dot(left_normal);

			T b11 = a11 * a11 * Lb + 2.0 * a11 * a21 * Mb + a21 * a21 * Nb - La;
			T b12 = a11 * a12 * Lb + 2.0 * (a11 * a22 + a21 * a12) * Mb + a21 * a22 * Nb - Ma;
			T b22 = a12 * a12 * Lb + 2.0 * (a12 * a22) * Mb + a22 * a22 * Nb - Na;

			bool is_reverse{ false };
			if (std::abs(b22) > std::abs(b11))
			{
				std::swap(b22, b11);
				is_reverse = true;
			}

			T delta = b11 * b22 - b12 * b12;

			auto check_result = [&](const T real_root) -> bool
				{
					T result = real_root * real_root * b11 * b11 + 2.0 * b12 * real_root + b22;
					if (is_equal(result, 0.0, b11 * b11 * PRECISION<T>::value * 1e-3))
					{
						return true;
					}
					return false;
				};

			std::vector<T> real_roots = solve_plolynomial_2_degree<T>::solve(b11, b12, b22);
			if (real_roots.size() == 2)
			{
				// check middle point
				T test_value = (real_roots[0] + real_roots[1]) / 2.0;
				if (false == check_result(test_value))
				{
					return ENUM_NURBS::NURBS_ERROR;
				}
				real_roots.clear();
				real_roots.push_back(test_value);
			}

			if (is_reverse == false)
			{
				for (auto& root : real_roots)
				{
					Eigen::Vector<T, dim> tangent = root * left_ders(1, 0) + right_ders(1, 0);
					if (tangent.squaredNorm() < PRECISION<T>::value * 1e-6)
					{
						return ENUM_NURBS::NURBS_DEGENERATE;
					}
					tangent.normalize();
					tangents.push_back(tangent);
				}
			}
			else
			{
				for (auto& root : real_roots)
				{
					Eigen::Vector<T, dim> tangent = left_ders(1, 0) + root * right_ders(1, 0);
					if (tangent.squaredNorm() < PRECISION<T>::value * 1e-6)
					{
						return ENUM_NURBS::NURBS_DEGENERATE;
					}
					tangent.normalize();
					tangents.push_back(tangent);
				}
			}

			// check result
			for (auto& tangent : tangents)
			{
				T left_dot = left_normal.dot(tangent);
				if (left_dot > PRECISION<T>::value * 1e2)
				{
					return ENUM_NURBS::NURBS_ERROR;
				}
				T right_dot = right_normal.dot(tangent);
				if (right_dot > PRECISION<T>::value * 1e2)
				{
					return ENUM_NURBS::NURBS_ERROR;
				}
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		}

		static ENUM_NURBS evaluate_curvature_normal(T u, T v, T s, T t, const surf_type& left_surf, const surf_type& right_surf, Eigen::Vector<T, dim>& kvec)
		{
			return ENUM_NURBS::NURBS_SUCCESS;
		}
	};


	template<typename Scalar, int rows, int cols>
	struct Newton_Iter
	{
		virtual ENUM_NURBS evaluate(const Eigen::Vector<Scalar, cols>& init_param, Eigen::Vector<Scalar, cols>& function_value, Eigen::Matrix<Scalar, rows, cols>& Jacobi_mat) = 0;

		virtual ENUM_NURBS update_param(const Eigen::Vector<Scalar, cols>& current_param, const Eigen::Vector<Scalar, cols>& delta_vec, Eigen::Vector<Scalar, cols>& next_param) = 0;

		ENUM_NURBS solve(const Eigen::Vector<Scalar, cols>& init_param, Eigen::Vector<Scalar, cols>& result_param, const Scalar tol, const int max_iter_num)
		{
			Eigen::Vector<Scalar, rows> vec;
			Scalar min_distance = 100000;
			Eigen::Vector<Scalar, cols> current_param = init_param;
			Scalar tol2 = tol * tol;
			Eigen::Matrix<Scalar, rows, cols> mat;
			ENUM_NURBS err_code = ENUM_NURBS::NURBS_SUCCESS;
			//����������Ҫ�޸�
			for (int loop_index = 0; loop_index < max_iter_num; ++loop_index)
			{
				err_code = evaluate(current_param, vec, mat);
				if (err_code != ENUM_NURBS::NURBS_SUCCESS)
				{
					return err_code;
				}
				Scalar distance = vec.squaredNorm();
				if (distance < min_distance)
				{
					min_distance = distance;
					result_param = current_param;
					if (distance < tol2)
						return ENUM_NURBS::NURBS_SUCCESS;
				}


				Eigen::JacobiSVD<Eigen::Matrix<Scalar, rows, cols>, Eigen::ComputeThinU | Eigen::ComputeThinV> matSvd(mat);
				Eigen::Vector<Scalar, cols> delta = matSvd.solve(vec);
				if (matSvd.info() != Eigen::Success)
					return ENUM_NURBS::NURBS_ERROR;

				Eigen::Vector<Scalar, cols> next_param;
				err_code = update_param(current_param, -delta, next_param);
				if (err_code != ENUM_NURBS::NURBS_SUCCESS)
				{
					return err_code;
				}
				current_param = next_param;
			}
			return ENUM_NURBS::NURBS_ERROR;
		}

		virtual ~Newton_Iter() {};
	};

	template<typename Scalar, int dim>
	struct plane_2surf_iter : public Newton_Iter<Scalar, 2 * dim, 6>
	{
		static constexpr int param_count = 6;
		static constexpr int rows = 2 * dim;
		
		using surf_type = surface<Scalar, dim>;
		using vec_type = Eigen::Vector<Scalar, dim>;
		using param_type = Eigen::Vector<Scalar, param_count>;
		using value_type = Eigen::Vector<Scalar, rows>;
		using matrix_type = Eigen::Matrix<Scalar, rows, param_count>;

	private:
		surf_type* m_surf1;
		surf_type* m_surf2;

		vec_type m_origin;
		vec_type m_u_dir;
		vec_type m_v_dir;

		//ǰ����surf1 domain, ��������surf2 domain
		Box<Scalar, 4> m_surfs_domain;

		//ǰ����surf1��uv������, ��������surf2 uv������
		std::array<bool, 4> m_surfs_periods;
	public:

		plane_2surf_iter(surf_type* surf1, surf_type* surf2, const vec_type& origin, const vec_type& u_dir, const vec_type& v_dir) :
			m_surf1(surf1), m_surf2(surf2), m_origin(origin), m_u_dir(u_dir), m_v_dir(v_dir)
		{
			Box<Scalar, 2> surf_domain = surf1->get_domain();
			m_surfs_domain.Min.head<2>() = surf_domain.Min;
			m_surfs_domain.Max.head<2>() = surf_domain.Max;
			surf_domain = surf2->get_domain();
			m_surfs_domain.Min.tail<2>() = surf_domain.Min;
			m_surfs_domain.Max.tail<2>() = surf_domain.Max;
			
			std::array<bool, 2> surf_periods;
			
			surf1->is_period(surf_periods);
			std::copy(surf_periods.begin(), surf_periods.end(), m_surfs_periods.begin());

			surf2->is_period(surf_periods);
			std::copy(surf_periods.begin(), surf_periods.end(), m_surfs_periods.begin() + 2);

		}

		// (init_param[0], init_param[1]) = S1(u, v), (init_param[2], init_param[3]) = S2(u, v), (init_param[4], init_param[5]) = plane(u, v)
		// function_value[0,..., dim - 1] = surf1 - plane; function[dim, 2 * dim - 1] = surf2 - plane
		virtual ENUM_NURBS evaluate(const param_type& init_param, value_type& function_value, matrix_type& Jacobi_mat) override
		{
			// Eigen::Matrix<vec_type, 2, 2> ders1, ders2;
			// ENUM_NURBS err_code = m_surf1->derivative_on_surface<1>(init_param[0], init_param[1], ders1);
			
			Eigen::MatrixX<vec_type> ders1, ders2;
			ENUM_NURBS err_code = m_surf1->derivative_on_surface(1, init_param[0], init_param[1], ders1);
			if (err_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				return err_code;
			}
			// err_code = m_surf2->derivative_on_surface<1>(init_param[2], init_param[3], ders2);
			err_code = m_surf2->derivative_on_surface(1, init_param[2], init_param[3], ders2);
			if (err_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				return err_code;
			}

			vec_type point3 = init_param[4] * m_u_dir + init_param[5] * m_v_dir + m_origin;

			function_value.head<dim>() = ders1(0, 0) - point3;
			function_value.tail<dim>() = ders2(0, 0) - point3;
			Jacobi_mat.setZero();

			Jacobi_mat.template block<dim, 1>(0, 0) = ders1(1, 0);
			Jacobi_mat.template block<dim, 1>(0, 1) = ders1(0, 1);
			Jacobi_mat.template block<dim, 1>(0, 4) = -m_u_dir;
			Jacobi_mat.template block<dim, 1>(0, 5) = -m_v_dir;
			Jacobi_mat.template block<dim, 1>(dim, 2) = ders2(1, 0);
			Jacobi_mat.template block<dim, 1>(dim, 3) = ders2(0, 1);
			Jacobi_mat.template block<dim, 1>(dim, 4) = -m_u_dir;
			Jacobi_mat.template block<dim, 1>(dim, 5) = -m_v_dir;
			return ENUM_NURBS::NURBS_SUCCESS;

		}
		virtual ENUM_NURBS update_param(const param_type& current_param, const param_type& delta_vec, param_type& next_param) override
		{
			auto modify_param = [](const Scalar low, const Scalar high, Scalar& param)
				{
					if (param < low)
					{
						param = low;
					}
					else if (param > high)
					{
						param = high;
					}
				};
			
			// ��param��һ��һ����׼������
			auto regular_param = [](const Scalar low, const Scalar high, Scalar& param)
				{
					Scalar period = high - low;
					if (param > high)
					{
						int k = static_cast<int>((param - low) / period);
						param -= k * period;
					}
					else if (param < low)
					{
						int k = static_cast<int>((high - param) / period);
						param += k * period;
					}
					return;
				};
			next_param = current_param + delta_vec;
			Scalar ratio = 1.0;
			bool need_scale{ false };
			for (size_t index = 0; index < 4; ++index)
			{
				Scalar low = m_surfs_domain.Min[index];
				Scalar high = m_surfs_domain.Max[index];
				if (compare<Scalar>::is_less(next_param[index], low, TDEFAULT_ERROR<Scalar>::value)
					|| compare<Scalar>::is_greater(next_param[index], high, TDEFAULT_ERROR<Scalar>::value))
				{
					if (m_surfs_periods[index] == true)
					{
						regular_param(low, high, next_param[index]);
					}
					else
					{
						need_scale = true;
						if (delta_vec[index] > 0)
						{
							Scalar len = high - current_param[index];
							len /= delta_vec[index];
							ratio = std::min(len, ratio);
						}
						else
						{
							Scalar len = low - current_param[index];
							len /= delta_vec[index];
							ratio = std::min(len, ratio);
						}
					}
				}
			}
		
			// assert(ratio > 0);
			// if (need_scale == true)
			// {
			// 	next_param = current_param + ratio * delta_vec;
			// }
			for (int index = 0; index < 4; ++index)
			{
				modify_param(m_surfs_domain.Min[index], m_surfs_domain.Max[index], next_param[index]);
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		}
	};


	template<typename Scalar, int dim>
	class int_curve : public curve<Scalar, dim>
	{
		static_assert(dim == 3, "dim != 3");
		
		using surf_type = surface<Scalar, dim>;
		
	private:
		surfs_int_points_chat<Scalar, dim> m_int_chat;
		std::vector<Scalar> m_params;
		// std::vector<Eigen::Vector<Scalar, dim>> m_ders;

		std::vector<Scalar> m_ms;
		std::vector<Scalar> m_ns;
		std::vector<Scalar> m_ps;
		std::vector<Scalar> m_qs;
		
		surf_type* m_left_surf;
		surf_type* m_right_surf;


		int_curve() = default;
	public:

		static int_curve* create(const surfs_int_points_chat<Scalar, dim>& int_chat, surf_type* left_surf, surf_type* right_surf)
		{
			int_curve* icurve = new int_curve();
			
			icurve->m_int_chat = int_chat;
			icurve->m_left_surf = left_surf;
			icurve->m_right_surf = right_surf;
			size_t points_count = int_chat.m_inter_points.size();
			icurve->m_params.reserve(points_count);

			std::vector<Scalar> as;
			std::vector<Scalar> bs;

			as.reserve(points_count - 1);
			bs.reserve(points_count - 1);
			Scalar current_param{ 0.0 };
			icurve->m_params.push_back(current_param);

			if (points_count <= 1)
			{
				return nullptr;
			}
			as.push_back(1.0);
			Eigen::Vector<Scalar, dim> currnet_normal = int_chat.m_inter_points[1].m_point - int_chat.m_inter_points[0].m_point;
			Scalar normal_len = currnet_normal.norm();
			current_param += normal_len;
			icurve->m_params.push_back(current_param);
			currnet_normal /= normal_len;
			for (size_t index = 1; index < points_count - 1; ++index)
			{
				const auto& current_point = int_chat.m_inter_points[index].m_point;
				const auto& next_point = int_chat.m_inter_points[index + 1].m_point;
				const auto& current_tangent = int_chat.m_inter_points[index].m_tangent;
				
				Eigen::Vector<Scalar, dim> next_normal = next_point - current_point;
				normal_len = next_normal.norm();
				current_param += normal_len;
				icurve->m_params.push_back(current_param);
				next_normal /= normal_len;
				Scalar e = currnet_normal.dot(current_tangent) / next_normal.dot(current_tangent);
				if (std::isnan(e))
				{
					assert(false);
					return nullptr;
				}
				if (e > 1.0)
				{
					as.push_back(1.0 / e);
					bs.push_back(1.0);
				}
				else
				{
					as.push_back(1.0);
					bs.push_back(e);
				}
				currnet_normal = next_normal;
			}
	
			const auto& start_point = int_chat.m_inter_points.front().m_point;
			const auto& end_point = int_chat.m_inter_points.back().m_point;
			const auto diff_vec = start_point - end_point;

			if (diff_vec.squaredNorm() > PRECISION<Scalar>::value * 1e-6)
			{
				//open curve
				bs.push_back(1.0);
			}
			else
			{
				// closed curve
				const auto& start_tangent = int_chat.m_inter_points.front().m_tangent;
				const auto& end_tangent = int_chat.m_inter_points.back().m_tangent;
				
				const auto& second_point = int_chat.m_inter_points[1].m_point;
				const auto& seond_tangent = int_chat.m_inter_points[1].m_tangent;

				const auto& third_point = int_chat.m_inter_points[points_count - 2].m_point;
				const auto& third_tangent = int_chat.m_inter_points[points_count - 2].m_tangent;

				Eigen::Vector<Scalar, dim> first_normal = second_point - start_point;
				Eigen::Vector<Scalar, dim> last_normal = end_point - third_point;
				Scalar en = first_normal.dot(start_tangent) / last_normal.dot(end_tangent);
				if (std::isnan(en))
				{
					assert(false);
					return nullptr;
				}
				if (en > 1.0)
				{
					as[0] = 1.0 / en;
					bs.push_back(1.0);
				}
				else
				{
					as[0] = 1.0;
					bs.push_back(en);
				}
			}

		
			icurve->m_ms.reserve(points_count - 1);
			icurve->m_ns.reserve(points_count - 1);
			icurve->m_ps.reserve(points_count - 1);
			icurve->m_qs.reserve(points_count - 1);
			for (size_t index = 0; index < points_count - 1; ++index)
			{
				icurve->m_qs.push_back(icurve->m_params[index]);
				icurve->m_ps.push_back(as[index]);
				Scalar d = icurve->m_params[index + 1] - icurve->m_params[index];
				icurve->m_ns.push_back((3.0 - 2.0 * as[index] - bs[index]) / d);
				icurve->m_ms.push_back((bs[index] + as[index] - 2.0) / (d * d));
			}
			return icurve;
		}

		
		Interval<Scalar> get_interval() const
		{
			Interval<Scalar> interval(m_params.front(), m_params.back());
			return interval;
		}

		ENGEOMETRYTYPE get_type() const { return ENGEOMETRYTYPE::NURBS_INT_CURVE; }

		virtual ENUM_NURBS point_on_curve(Scalar u, Eigen::Vector<Scalar, dim>& point) const override
		{
			size_t param_count = m_params.size();
			if (param_count < 2)
			{
				return ENUM_NURBS::NURBS_ERROR;
			}
			if (u < m_params.front() || u > m_params.back())
			{
				return ENUM_NURBS::NURBS_ERROR;
			}

			size_t index = 0;
			for (; index < param_count - 1; ++index)
			{
				if (m_params[index] <= u && u <= m_params[index + 1])
				{
					break;
				}
			}

			Scalar delta = u - m_params[index];
			Scalar delta2 = delta * delta;
			Scalar param = m_ms[index] * delta2 * delta + m_ns[index] * delta2 + m_ps[index] * delta + m_qs[index];
			ENUM_NURBS err_code = point_on_segment(param, index, point);
			if (err_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				return err_code;
			}
			return ENUM_NURBS::NURBS_SUCCESS;
		};
	private:

		// m_params[start_index] <= u <= m_params[start_index + 1]
		ENUM_NURBS point_on_segment(const Scalar u, const size_t start_index, Eigen::Vector<Scalar, dim>& point) const
		{
			Scalar u_vec = u - m_params[start_index];
			Scalar param_vec = m_params[start_index + 1] - m_params[start_index];
			const auto& start_point = m_int_chat.m_inter_points[start_index].m_point;
			const auto& next_point = m_int_chat.m_inter_points[start_index + 1].m_point;
			
			const auto& start_param = m_int_chat.m_inter_points[start_index].m_uv;
			const auto& next_param = m_int_chat.m_inter_points[start_index + 1].m_uv;
			
			Eigen::Vector<Scalar, dim> diff_vec = next_point - start_point;
			Scalar len_ratio = u_vec / param_vec;
			if (std::isnan(len_ratio))
			{
				assert(false);
				return ENUM_NURBS::NURBS_ERROR;
			}
			Eigen::Vector<Scalar, dim> origin = len_ratio * diff_vec + start_point;
			const Eigen::Vector<Scalar, 4> estimat_param = len_ratio * (next_param - start_param) + start_param;
			std::array<Eigen::Vector<Scalar, dim>, 2> lcs;
			diff_vec.normalize();
			vec3_utils<Scalar>::create_lcs(diff_vec, lcs);
			plane_2surf_iter<Scalar, dim> iter_func(m_left_surf, m_right_surf, origin, lcs[0], lcs[1]);
			Eigen::Vector<Scalar, 6> init_param{estimat_param[0], estimat_param[1], estimat_param[2], estimat_param[3], 0, 0};
			Eigen::Vector<Scalar, 6> result_param;;
			ENUM_NURBS err_code = iter_func.solve(init_param, result_param, PRECISION<Scalar>::value, 20);
			if (err_code != ENUM_NURBS::NURBS_SUCCESS)
			{
				return err_code;
			}
			m_left_surf->point_on_surface(result_param[0], result_param[1], point);
			return ENUM_NURBS::NURBS_SUCCESS;
		};
	};
}
