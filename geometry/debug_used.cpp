#include "debug_used.h"
#include <json.hpp>
#include <fstream>

namespace tnurbs
{
    void printEigenVector(const Eigen::VectorX<double> &vec)
    {
        std::cout << vec << std::endl;
    }
    void save_box(const std::vector<Box<double, 2>> &boxes, const char *path)
    {
        int cols = boxes.size();
        std::ofstream outfile2(path);
        outfile2.precision(16);

        for (int j = 0; j < cols; ++j)
        {
            Box<double, 2> box = boxes[j];
            
            outfile2 << "v " << box.Min[0] << " " <<
            box.Min[1] << " " << 0 << std::endl;
            outfile2 << "v " << box.Max[0] << " " <<
                box.Min[1] << " " << 0 << std::endl;
            outfile2 << "v " << box.Max[0] << " " <<
                box.Max[1] << " " << 0 << std::endl;
            outfile2 << "v " << box.Min[0] << " " <<
                box.Max[1] << " " << 0 << std::endl;
        }
        for (int j = 0; j < cols; ++j)
        {
            outfile2 << "l ";
            outfile2 << j * 4 + 1 << " " << j * 4 + 2 << std::endl;
            outfile2 << "l ";
            outfile2 << j * 4 + 2 << " " << j * 4 + 3 << std::endl;
            outfile2 << "l ";
            outfile2 << j * 4 + 3 << " " << j * 4 + 4 << std::endl;
            outfile2 << "l ";
            outfile2 << j * 4 + 4 << " " << j * 4 + 1 << std::endl;
        }
        outfile2.close();
    }

    void save_obj(Eigen::Matrix<double, 3, Eigen::Dynamic>& mat, const char* path)
    {
        int cols = mat.cols();
        std::ofstream outfile2(path);
        outfile2.precision(16);

        for (int j = 0; j < cols; ++j)
        {
            Eigen::Vector3d point = mat.col(j);

            outfile2 << "v " << point[0] << " " <<
                point[1] << " " << point[2] << std::endl;
        }
        outfile2.close();
    }
    
    void save_obj(std::vector<Eigen::Vector3d>& points, const char* path)
    {
        int cols = points.size();
        std::ofstream outfile2(path);
        outfile2.precision(16);

        for (int j = 0; j < cols; ++j)
        {
            Eigen::Vector3d point = points[j];

            outfile2 << "v " << point[0] << " " <<
                point[1] << " " << point[2] << std::endl;
        }
        outfile2.close();
    }
    
    void save_obj(const Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> &mat, const char *path)
    {
        int rows = mat.rows();
        int cols = mat[0].cols();
        std::ofstream outfile2(path);

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                Eigen::Vector3d point = mat[i].col(j);
                outfile2 << "v " << point[0] << " " <<
                point[1] << " " << point[2] << std::endl;
            }
        }
        outfile2.close();
    }

    void save_obj(Eigen::Matrix<double, 2, Eigen::Dynamic>& mat, const char* path)
    {
        int cols = mat.cols();
        std::ofstream outfile2(path);
        outfile2.precision(16);

        for (int j = 0; j < cols; ++j)
        {
            Eigen::Vector2d point = mat.col(j);

            outfile2 << "v " << point[0] << " " <<
                point[1] << " " << 0 << std::endl;
        }
        outfile2.close();
    }

    void printEigenMatrix(const Eigen::Matrix<double, 2, Eigen::Dynamic> &mat)
    {
        std::cout << mat << std::endl;
    }

    void printEigenMatrix(const Eigen::Matrix<double, 3, Eigen::Dynamic> &mat)
    {
        std::cout << mat << std::endl;
    }
    void printEigenMatrix(const Eigen::SparseMatrix<double> &mat)
    {
        std::cout << mat << std::endl;
    }

    void printEigenMatrix(const Eigen::VectorX<Eigen::Matrix<double, 3, Eigen::Dynamic>> &mat)
    {
        int cols = mat.rows();
        for (int index = 0; index < cols; ++index)
        {
            std::cout << "********" << std::endl;
            std::cout << mat[index] << std::endl;
            std::cout << "********" << std::endl;
        }
        
    }

    void printEigenMatrix(const Eigen::Matrix<double, 4, Eigen::Dynamic> &mat)
    {
        std::cout << mat << std::endl;
    }

    void printEigenMatrix(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &mat)
    {
        std::cout << mat << std::endl;
    }

    // template<typename T, int dim, bool is_rational>
    ENUM_NURBS save_obj(const nurbs_curve<double, 2, false, -1, -1> &nurbs_cur, const char *path)
    {
        std::array<double, 2> end_konts;
        nurbs_cur.get_ends_knots(end_konts); 
        double step = (end_konts[1] - end_konts[0]) / 100.0;
        std::vector<Eigen::Vector<double, 2>> points;
        for (int u_index = 0; u_index < 100; ++u_index)
        {
            Eigen::Vector<double, 2> point;
            nurbs_cur.point_on_curve(u_index * step + end_konts[0], point);
            points.push_back(point);
        }
        std::ofstream outfile2(path);
        for (auto point : points)
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << 0 << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_obj(const nurbs_curve<double, 2, true, -1, -1> &nurbs_cur, const char *path)
    {
        std::array<double, 2> end_konts;
        nurbs_cur.get_ends_knots(end_konts); 
        double step = (end_konts[1] - end_konts[0]) / 100.0;
        std::vector<Eigen::Vector<double, 2>> points;
        for (int u_index = 0; u_index < 100; ++u_index)
        {
            Eigen::Vector<double, 2> point;
            nurbs_cur.point_on_curve(u_index * step + end_konts[0], point);
            points.push_back(point);
        }
        std::ofstream outfile2(path);
        for (auto point : points)
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << 0 << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_obj(const nurbs_curve<double, 3, false, -1, -1> &nurbs_cur, const char *path)
    {
        std::array<double, 2> end_konts;
        nurbs_cur.get_ends_knots(end_konts); 
        double step = (end_konts[1] - end_konts[0]) / 100.0;
        std::vector<Eigen::Vector<double, 3>> points;
        for (int u_index = 0; u_index < 100; ++u_index)
        {
            Eigen::Vector<double, 3> point;
            nurbs_cur.point_on_curve(u_index * step + end_konts[0], point);
            points.push_back(point);
        }
        std::ofstream outfile2(path);
        for (auto point : points)
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << point[2] << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_obj(const nurbs_curve<double, 3, true, -1, -1> &nurbs_cur, const char *path)
    {
        std::array<double, 2> end_konts;
        nurbs_cur.get_ends_knots(end_konts); 
        double step = (end_konts[1] - end_konts[0]) / 100.0;
        std::vector<Eigen::Vector<double, 3>> points;
        for (int u_index = 0; u_index < 100; ++u_index)
        {
            Eigen::Vector<double, 3> point;
            nurbs_cur.point_on_curve(u_index * step + end_konts[0], point);
            points.push_back(point);
        }
        std::ofstream outfile2(path);
        for (auto point : points)
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << point[2] << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }


    ENUM_NURBS save_obj(const nurbs_surface<double, 2, -1, -1, -1, -1, false> &surf, const char *path)
    {
        Eigen::VectorX<double> u_knots_vector = surf.get_u_knots();
        Eigen::VectorX<double> v_knots_vector = surf.get_v_knots();;
        int u_knots_size = u_knots_vector.size();
        int v_knots_size = v_knots_vector.size();
        double u_high = u_knots_vector[u_knots_size - 1];
        double u_low = u_knots_vector[0];
        double v_high = v_knots_vector[v_knots_size - 1];
        double v_low = v_knots_vector[0];
        double u_step = (u_high - u_low) / 100.0;
        double v_step = (v_high - v_low) / 100.0;
        std::vector<Eigen::Vector<double, 2>> points;
        for (int u_index = 0; u_index < 100; ++u_index)
        {
            for (int v_index = 0; v_index < 100; ++v_index)
            {
                Eigen::Vector<double, 2> point;
                surf.point_on_surface(u_index * u_step + u_low, v_index * v_step + v_low, point);
                points.push_back(point);
            }
        }
        std::ofstream outfile2(path);
        for (auto point : points)
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << point[2] << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_obj(const nurbs_surface<double, 2, -1, -1, -1, -1, true> &surf, const char *path)
    {
        Eigen::VectorX<double> u_knots_vector = surf.get_u_knots();
        Eigen::VectorX<double> v_knots_vector = surf.get_v_knots();;
        int u_knots_size = u_knots_vector.size();
        int v_knots_size = v_knots_vector.size();
        double u_high = u_knots_vector[u_knots_size - 1];
        double u_low = u_knots_vector[0];
        double v_high = v_knots_vector[v_knots_size - 1];
        double v_low = v_knots_vector[0];
        double u_step = (u_high - u_low) / 100.0;
        double v_step = (v_high - v_low) / 100.0;
        std::vector<Eigen::Vector<double, 2>> points;
        for (int u_index = 0; u_index < 100; ++u_index)
        {
            for (int v_index = 0; v_index < 100; ++v_index)
            {
                Eigen::Vector<double, 2> point;
                surf.point_on_surface(u_index * u_step + u_low, v_index * v_step + v_low, point);
                points.push_back(point);
            }
        }
        std::ofstream outfile2(path);
        for (auto point : points)
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << point[2] << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_obj(const nurbs_surface<double, 3, -1, -1, -1, -1, false> &surf, const char *path)
    {
        Eigen::VectorX<double> u_knots_vector = surf.get_u_knots();
        Eigen::VectorX<double> v_knots_vector = surf.get_v_knots();;
        int u_knots_size = u_knots_vector.size();
        int v_knots_size = v_knots_vector.size();
        double u_high = u_knots_vector[u_knots_size - 1];
        double u_low = u_knots_vector[0];
        double v_high = v_knots_vector[v_knots_size - 1];
        double v_low = v_knots_vector[0];
        double u_step = (u_high - u_low) / 100.0;
        double v_step = (v_high - v_low) / 100.0;
        std::vector<Eigen::Vector<double, 3>> points;
        for (int u_index = 0; u_index < 100; ++u_index)
        {
            for (int v_index = 0; v_index < 100; ++v_index)
            {
                Eigen::Vector<double, 3> point;
                surf.point_on_surface(u_index * u_step + u_low, v_index * v_step + v_low, point);
                points.push_back(point);
            }
        }
        std::ofstream outfile2(path);
        for (auto point : points)
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << point[2] << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_obj(const nurbs_surface<double, 3, -1, -1, -1, -1, true> &surf, const char *path)
    {
        Eigen::VectorX<double> u_knots_vector = surf.get_u_knots();
        Eigen::VectorX<double> v_knots_vector = surf.get_v_knots();;
        int u_knots_size = u_knots_vector.size();
        int v_knots_size = v_knots_vector.size();
        double u_high = u_knots_vector[u_knots_size - 1];
        double u_low = u_knots_vector[0];
        double v_high = v_knots_vector[v_knots_size - 1];
        double v_low = v_knots_vector[0];
        double u_step = (u_high - u_low) / 100.0;
        double v_step = (v_high - v_low) / 100.0;
        std::vector<Eigen::Vector<double, 3>> points;
        for (int u_index = 0; u_index < 100; ++u_index)
        {
            for (int v_index = 0; v_index < 100; ++v_index)
            {
                Eigen::Vector<double, 3> point;
                surf.point_on_surface(u_index * u_step + u_low, v_index * v_step + v_low, point);
                points.push_back(point);
            }
        }
        std::ofstream outfile2(path);
        for (auto point : points)
        {
            outfile2 << "v " << point[0] << " " <<
            point[1] << " " << point[2] << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_obj(const mesh<2> *surf, const char *path)
    {
        std::ofstream outfile2(path);
        for (auto point : surf->m_ders)
        {
            outfile2 << "v " << point(0, 0)[0] << " " <<
            point(0, 0)[1] << " " << point(0, 0)[2] << std::endl;
        }
        for (auto index : surf->m_indexs)
        {
            int count = index.size();
            outfile2 << "f ";
            for (int i = 0; i < count; ++i)
            {
                outfile2 << index[i] + 1 << " ";
            }
            outfile2 << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    ENUM_NURBS save_obj(const mesh<1>* cur, const char* path)
    {
        std::ofstream outfile2(path);
        for (auto point : cur->m_ders)
        {
            outfile2 << "v " << point(0, 0)[0] << " " <<
                point(0, 0)[1] << " " << point(0, 0)[2] << std::endl;
        }
        for (auto index : cur->m_indexs)
        {
            int count = index.size();
            outfile2 << "l ";
            for (int i = 0; i < count; ++i)
            {
                outfile2 << index[i] + 1 << " ";
            }
            outfile2 << std::endl;
        }
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_obj2(const nurbs_surface<double, 3, -1, -1, -1, -1, false>& surf, const char* path)
    {
        surface_mesh_helper<nurbs_surface<double, 3, -1, -1, -1, -1, false>> mh;
        disc_surface(&surf, mh, TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
        mesh<2> surface_mesh;
        mesh_help_to_mesh(mh, surface_mesh);

        save_obj(&surface_mesh, path);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_obj2(const nurbs_surface<double, 3, -1, -1, -1, -1, true>& surf, const char* path)
    {
        surface_mesh_helper<nurbs_surface<double, 3, -1, -1, -1, -1, true>> mh;
        disc_surface(&surf, mh, TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
        mesh<2> surface_mesh;
        mesh_help_to_mesh(mh, surface_mesh);

        save_obj(&surface_mesh, path);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    ENUM_NURBS save_obj2(const nurbs_curve<double, 3, false, -1, -1>& curv, const char* path)
    {
        curve_mesh_helper<nurbs_curve<double, 3, false, -1, -1>> mh;
        disc_curve(&curv, mh, TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
        mesh<1> curve_mesh;
        mesh_help_to_mesh(mh, curve_mesh);

        save_obj(&curve_mesh, path);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_obj2(const nurbs_curve<double, 3, true, -1, -1>& curv, const char* path)
    {
        curve_mesh_helper<nurbs_curve<double, 3, true, -1, -1>> mh;
        disc_curve(&curv, mh, TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
        mesh<1> curve_mesh;
        mesh_help_to_mesh(mh, curve_mesh);

        save_obj(&curve_mesh, path);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    ENUM_NURBS save_chat_points(const std::vector<surfs_int_points_chat<double, 3>>& points_chats, const char* path)
    {
        nlohmann::json chat_json;
        chat_json["Name"] = "surface_int_point_chats";
        chat_json["CurveCount"] = points_chats.size();
        int index = 0;
        for (const auto& chat : points_chats)
        {
            chat_json["Curve" + std::to_string(index)] = {};
            nlohmann::json& curve_data = chat_json["Curve" + std::to_string(index)];
            curve_data["PointsCount"] = chat.m_inter_points.size();
			curve_data["Point3d"] = {};
			curve_data["Param"] = {};
            curve_data["Transversal"] = chat.m_is_transversal;
            for (const auto& int_point : chat.m_inter_points)
            {
				curve_data["Point3d"].push_back(int_point.m_point[0]);
				curve_data["Point3d"].push_back(int_point.m_point[1]);
				curve_data["Point3d"].push_back(int_point.m_point[2]);

				curve_data["Param"].push_back(int_point.m_uv[0]);
				curve_data["Param"].push_back(int_point.m_uv[1]);
				curve_data["Param"].push_back(int_point.m_uv[2]);
				curve_data["Param"].push_back(int_point.m_uv[3]);
            }
            index += 1;
        }
        std::ofstream outfile2(path);
        outfile2 << chat_json;
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS read_chat_points(std::vector<surfs_int_points_chat<double, 3>>& points_chat, const char* path)
    {
        nlohmann::json chat_json;
        std::ifstream jfile(path);
        jfile >> chat_json;
        size_t curve_count = chat_json.at("CurveCount").get<size_t>();

        points_chat.resize(curve_count);
        for (size_t index = 0; index < curve_count; ++index)
        {
            surfs_int_points_chat<double, 3>& chat = points_chat[index];
            nlohmann::json curve_data = chat_json["Curve" + std::to_string(index)];
            nlohmann::json points_data = curve_data["Point3d"];
            nlohmann::json param_data = curve_data["Param"];
            size_t points_count = curve_data["PointsCount"].get<size_t>();
            chat.m_inter_points.resize(points_count);
            chat.m_is_transversal = curve_data["Transversal"].get<bool>();
            for (size_t i = 0; i < points_count; ++i)
            {
                size_t k = i * 3;
                chat.m_inter_points[i].m_point = Eigen::Vector3d(points_data[k].get<double>(), points_data[k + 1].get<double>(), points_data[k + 2].get<double>());
                size_t j = i * 4;
                chat.m_inter_points[i].m_uv = Eigen::Vector4d(param_data[j].get<double>(), param_data[1 + j].get<double>(), param_data[2 + j].get<double>(), param_data[3 + j].get<double>());
            }
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    ENUM_NURBS read_chat_points(surf_surf_int<double, 3>& int_data, const std::string& path)
    {
        nlohmann::json chat_json;
        std::ifstream jfile(path);
        jfile >> chat_json;
        size_t curve_count = chat_json.at("CurveCount").get<size_t>();
        std::vector<surfs_int_points_chat<double, 3>>& points_chat = int_data.m_int_chats;
        points_chat.resize(curve_count);
        for (size_t index = 0; index < curve_count; ++index)
        {
            surfs_int_points_chat<double, 3>& chat = points_chat[index];
            nlohmann::json curve_data = chat_json["Curve" + std::to_string(index)];
            nlohmann::json points_data = curve_data["Point3d"];
            nlohmann::json param_data = curve_data["Param"];
            size_t points_count = curve_data["PointsCount"].get<size_t>();
            chat.m_inter_points.resize(points_count);
            chat.m_is_transversal = curve_data["Transversal"].get<bool>();
            for (size_t i = 0; i < points_count; ++i)
            {
                size_t k = i * 3;
                chat.m_inter_points[i].m_point = Eigen::Vector3d(points_data[k].get<double>(), points_data[k + 1].get<double>(), points_data[k + 2].get<double>());
                size_t j = i * 4;
                chat.m_inter_points[i].m_uv = Eigen::Vector4d(param_data[j].get<double>(), param_data[1 + j].get<double>(), param_data[2 + j].get<double>(), param_data[3 + j].get<double>());
            }
        }
        size_t isolate_point_count = chat_json["IsolatePointCount"].get<size_t>();
        int_data.m_isolate_points.resize(isolate_point_count);
        nlohmann::json isolate_points = chat_json["IsolatePoint"];
        nlohmann::json isolate_param = chat_json["IsolateParam"];
        for (size_t index = 0; index < isolate_point_count; ++index)
        {
            size_t k = index * 3;
            int_data.m_isolate_points[index].m_point = Eigen::Vector3d(isolate_points[k].get<double>(), isolate_points[k + 1].get<double>(), isolate_points[k + 2].get<double>());
            size_t j = index * 4;
            int_data.m_isolate_points[index].m_uv = Eigen::Vector4d(isolate_param[j].get<double>(), isolate_param[j + 1].get<double>(), isolate_param[j + 2].get<double>(), isolate_param[j + 3].get<double>());
        }
        return ENUM_NURBS::NURBS_SUCCESS;
    }

    ENUM_NURBS save_chat_points(const surf_surf_int<double, 3>& points_chat, const char* path)
    {
        nlohmann::json chat_json;
        const std::vector<surfs_int_points_chat<double, 3>>& curves_data = points_chat.m_int_chats;
        chat_json["Name"] = "surface_int_point_chats";
        chat_json["CurveCount"] = curves_data.size();
        int index = 0;
        for (const auto& chat : curves_data)
        {
            chat_json["Curve" + std::to_string(index)] = {};
            nlohmann::json& curve_data = chat_json["Curve" + std::to_string(index)];
            curve_data["PointsCount"] = chat.m_inter_points.size();
			curve_data["Point3d"] = {};
			curve_data["Param"] = {};
            curve_data["Transversal"] = chat.m_is_transversal;
            for (const auto& int_point : chat.m_inter_points)
            {
				curve_data["Point3d"].push_back(int_point.m_point[0]);
				curve_data["Point3d"].push_back(int_point.m_point[1]);
				curve_data["Point3d"].push_back(int_point.m_point[2]);

				curve_data["Param"].push_back(int_point.m_uv[0]);
				curve_data["Param"].push_back(int_point.m_uv[1]);
				curve_data["Param"].push_back(int_point.m_uv[2]);
				curve_data["Param"].push_back(int_point.m_uv[3]);
            }
            index += 1;
        }
        const std::vector<surf_surf_intersect_point<double, 3>>& isolate_data = points_chat.m_isolate_points;
        chat_json["IsolateParam"] = {};
        chat_json["IsolatePoint"] = {};
        chat_json["IsolatePointCount"] = isolate_data.size();
        for (const auto& isolate : isolate_data)
        {
            chat_json["IsolatePoint"].push_back(isolate.m_point[0]);
            chat_json["IsolatePoint"].push_back(isolate.m_point[1]);
            chat_json["IsolatePoint"].push_back(isolate.m_point[2]);
            
            chat_json["IsolateParam"].push_back(isolate.m_uv[0]);
            chat_json["IsolateParam"].push_back(isolate.m_uv[1]);
            chat_json["IsolateParam"].push_back(isolate.m_uv[2]);
            chat_json["IsolateParam"].push_back(isolate.m_uv[3]);
        }
        std::ofstream outfile2(path);
        outfile2 << chat_json;
        outfile2.close();
        return ENUM_NURBS::NURBS_SUCCESS;
    }
    ENUM_NURBS save_chat_points_file(const surf_surf_int<double, 3>& points_chat, const char* path)
    {
        std::vector<Eigen::Vector3d> points;
        const std::vector<surfs_int_points_chat<double, 3>>& curves_data = points_chat.m_int_chats;
        for (const auto& chat : curves_data)
        {
            size_t points_count = chat.m_inter_points.size();
            for (const auto& int_point : chat.m_inter_points)
            {
                points.push_back(int_point.m_point);
            }
        }
        save_obj(points, path);
        return ENUM_NURBS::NURBS_SUCCESS;
    }

}


