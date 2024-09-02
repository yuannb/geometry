#include "debug_used.h"

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

    ENUM_NURBS save_obj2(const nurbs_curve<double, 3, false, -1, -1>& curv, const char* path)
    {
        curve_mesh_helper<nurbs_curve<double, 3, false, -1, -1>> mh;
        disc_curve(&curv, mh, TDEFAULT_ERROR<double>::value, 10.0, 0.1, 1.0);
        mesh<1> curve_mesh;
        mesh_help_to_mesh(mh, curve_mesh);

        save_obj(&curve_mesh, path);
        return ENUM_NURBS::NURBS_SUCCESS;
    }
}


