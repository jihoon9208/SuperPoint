#include <iostream>
#include <cstdio>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <numpy/ndarrayobject.h>
#include "boost/tuple/tuple.hpp"
#include "boost/python/object.hpp"
#include <boost/tuple/tuple_comparison.hpp>
#include <limits>
#include <map>
#include <cmath>
#include "connected_components.cpp"
#include "random_subgraph.cpp"
#include <fstream>
namespace bp = boost::python;
namespace ei = Eigen;
namespace bpn = boost::python::numpy;

typedef ei::Matrix<float, 3, 3> Matrix3f;
typedef ei::Matrix<float, 3, 1> Vector3f;
typedef ei::Matrix<double, 3, 1> Vector3d;
typedef ei::Matrix<double, 3, 3> Matrix3d;

typedef boost::tuple< std::vector< std::vector<float> >, std::vector< std::vector<uint8_t> >, std::vector< std::vector<uint32_t> >, std::vector<std::vector<uint32_t> > > Custom_tuple;
typedef boost::tuple< std::vector< std::vector<uint32_t> >, std::vector<uint32_t> > Components_tuple;
typedef boost::tuple< std::vector<uint8_t>, std::vector<uint8_t> > Subgraph_tuple;

typedef boost::tuple< uint32_t, uint32_t, uint32_t > Space_tuple;

struct VecToArray
{//converts a vector<uint8_t> to a numpy array
    static PyObject * convert(const std::vector<uint8_t> & vec) {
    npy_intp dims = vec.size();
    PyObject * obj = PyArray_SimpleNew(1, &dims, NPY_UINT8);
    void * arr_data = PyArray_DATA((PyArrayObject*)obj);
    memcpy(arr_data, &vec[0], dims * sizeof(uint8_t));
    return obj;
    }
};

template <class T>
struct VecvecToArray
{//converts a vector< vector<uint32_t> > to a numpy 2d array
    static PyObject * convert(const std::vector< std::vector<T> > & vecvec)
    {
        npy_intp dims[2];
        dims[0] = vecvec.size();
        dims[1] = vecvec[0].size();
        PyObject * obj;
        if (typeid(T) == typeid(uint8_t))
            obj = PyArray_SimpleNew(2, dims, NPY_UINT8);
        else if (typeid(T) == typeid(float))
            obj = PyArray_SimpleNew(2, dims, NPY_FLOAT32);
        else if (typeid(T) == typeid(uint32_t))
            obj = PyArray_SimpleNew(2, dims, NPY_UINT32);
        void * arr_data = PyArray_DATA((PyArrayObject*)obj);
        std::size_t cell_size = sizeof(T);
        for (std::size_t i = 0; i < dims[0]; i++)
        {
            memcpy(arr_data + i * dims[1] * cell_size, &(vecvec[i][0]), dims[1] * cell_size);
        }
        return obj;
    }
};

struct VecToArray32
{//converts a vector<uint32_t> to a numpy array
    static PyObject * convert(const std::vector<uint32_t> & vec)
    {
        npy_intp dims = vec.size();
        PyObject * obj = PyArray_SimpleNew(1, &dims, NPY_UINT32);
        void * arr_data = PyArray_DATA((PyArrayObject*)obj);
        memcpy(arr_data, &vec[0], dims * sizeof(uint32_t));
        return obj;
    }
};
template<class T>
struct VecvecToList
{//converts a vector< vector<T> > to a list
        static PyObject* convert(const std::vector< std::vector<T> > & vecvec)
    {
        boost::python::list* pylistlist = new boost::python::list();
        for(size_t i = 0; i < vecvec.size(); i++)
        {
            boost::python::list* pylist = new boost::python::list();
            for(size_t j = 0; j < vecvec[i].size(); j++)
            {
                pylist->append(vecvec[i][j]);
            }
            pylistlist->append((pylist, pylist[0]));
        }
        return pylistlist->ptr();
    }
};

struct to_py_tuple
{//converts to a python tuple
    static PyObject* convert(const Custom_tuple & c_tuple){
        bp::list values;

        PyObject * pyo1 = VecvecToArray<float>::convert(c_tuple.get<0>());
        PyObject * pyo2 = VecvecToArray<uint8_t>::convert(c_tuple.get<1>());
        PyObject * pyo3 = VecvecToArray<uint32_t>::convert(c_tuple.get<2>());
        PyObject * pyo4 = VecvecToArray<uint32_t>::convert(c_tuple.get<3>());

        values.append(bp::handle<>(bp::borrowed(pyo1)));
        values.append(bp::handle<>(bp::borrowed(pyo2)));
        values.append(bp::handle<>(bp::borrowed(pyo3)));
        values.append(bp::handle<>(bp::borrowed(pyo4)));

        return bp::incref( bp::tuple( values ).ptr() );
    }
};

struct to_py_tuple_components
{//converts output to a python tuple
    static PyObject* convert(const Components_tuple& c_tuple){
        bp::list values;
        //add all c_tuple items to "values" list

        PyObject * vecvec_pyo = VecvecToList<uint32_t>::convert(c_tuple.get<0>());
        PyObject * vec_pyo = VecToArray32::convert(c_tuple.get<1>());

        values.append(bp::handle<>(bp::borrowed(vecvec_pyo)));
        values.append(bp::handle<>(bp::borrowed(vec_pyo)));

        return bp::incref( bp::tuple( values ).ptr() );
    }
};

struct to_py_tuple_subgraph
{//converts output to a python tuple
    static PyObject* convert(const Subgraph_tuple& s_tuple){
        bp::list values;
        //add all c_tuple items to "values" list

        PyObject * vec_pyo1 = VecToArray::convert(s_tuple.get<0>());
        PyObject * vec_pyo2 = VecToArray::convert(s_tuple.get<1>());

        values.append(bp::handle<>(bp::borrowed(vec_pyo1)));
        values.append(bp::handle<>(bp::borrowed(vec_pyo2)));

        return bp::incref( bp::tuple( values ).ptr() );
    }
};

double four_normal(const float * xyz, const uint32_t * target, int i_ver, int nedg, int nornei, double xinorm, double weightk){
    
    ei::MatrixXf four_tposition(4,3);
    ei::Vector3d end_tposition;    
    std::size_t end_edg = nedg + 44;
    std::size_t end_nei = target[end_edg];

    end_tposition(0) = xyz[3 * end_nei];
    end_tposition(1) = xyz[3 * end_nei + 1];
    end_tposition(2) = xyz[3 * end_nei + 2];
     

    for (int i_nei = 0; i_nei < 4; i_nei++){
        nornei = target[nedg];
        four_tposition(i_nei,0) = xyz[3 * nornei];
        four_tposition(i_nei,1) = xyz[3 * nornei + 1];
        four_tposition(i_nei,2) = xyz[3 * nornei + 2];
        nedg++;
    }

    //compute four neighboor point normal estimation
    double four_norm = 0;
    for (std::size_t i_nei = 0; i_nei < 4; i_nei++) {
        double xj_norm = sqrtf(four_tposition(i_nei,0)*four_tposition(i_nei,0) + 
        four_tposition(i_nei,1)*four_tposition(i_nei,1) + four_tposition(i_nei,2)*four_tposition(i_nei,2));
        double weightj = exp(-0.2*fabsf(xj_norm-xinorm));
        ei::Vector3d e_position;     
        ei::Vector3d next_position;

        e_position(0) = weightk * (end_tposition(0) - xyz[3 * i_ver] );
        e_position(1) = weightk * (end_tposition(1) - xyz[3 * i_ver + 1] );
        e_position(2) = weightk * (end_tposition(2) - xyz[3 * i_ver + 2] );
        next_position(0) = weightj * (four_tposition(i_nei,0) - xyz[3 * i_ver] );
        next_position(1) = weightj * (four_tposition(i_nei,1) - xyz[3 * i_ver + 1] );
        next_position(2) = weightj * (four_tposition(i_nei,2) - xyz[3 * i_ver + 2] );
        
        four_norm += e_position.adjoint() * next_position;
        
    }
    
    return four_norm;
}

PyObject * estimate_geof(const bpn::ndarray & xyz ,const bpn::ndarray & target, int k_nn)
{//compute the following geometric features (geof) features of a point cloud:
 //linearity planarity scattering verticality
    std::size_t n_ver = bp::len(xyz);                                                     // xyz의 길이
    std::vector< std::vector< float > > geof(n_ver, std::vector< float >(7,0));             // geof 
    //--- read numpy array data---
    const uint32_t * target_data = reinterpret_cast<uint32_t*>(target.get_data());    
    const float * xyz_data = reinterpret_cast<float*>(xyz.get_data());
    
    std::size_t s_ver = 0;
    #pragma omp parallel for schedule(static)
    for (std::size_t i_ver = 0; i_ver < n_ver; i_ver++)
    {//each point can be treated in parallell independently
        //--- compute 3d covariance matrix of neighborhood ---
        ei::MatrixXf position(k_nn +1 ,3);
        ei::Vector3d end_position;    
        ei::MatrixXf four_position(4,3); 
                            
        std::size_t i_edg = k_nn * i_ver;
        std::size_t d_edg = k_nn * i_ver;
        std::size_t n_edg = k_nn * i_ver;
        std::size_t tn_edg = k_nn * i_ver;
        std::size_t ind_nei;
        std::size_t nor_nei;
        std::size_t tnor_nei;
        std::size_t end_edg = i_edg + 44;
        std::size_t end_nei = target_data[end_edg];

        end_position(0) = xyz_data[3 * end_nei];
        end_position(1) = xyz_data[3 * end_nei + 1];
        end_position(2) = xyz_data[3 * end_nei + 2];

        double xi_norm = sqrtf(xyz_data[3 * i_ver] * xyz_data[3 * i_ver] + 
            xyz_data[3 * i_ver+1] * xyz_data[3 * i_ver+1] + xyz_data[3 * i_ver+2] * xyz_data[3 * i_ver+2]);
        double xk_norm = sqrtf(end_position[0]*end_position[0] + 
            end_position[1]*end_position[1] + end_position[2]*end_position[2]);
        
        double weightk = exp(-0.2*abs(xk_norm-xi_norm));
        
        position(0,0) = xyz_data[3 * i_ver];
        position(0,1) = xyz_data[3 * i_ver + 1];
        position(0,2) = xyz_data[3 * i_ver + 2];

        for (std::size_t i_nei = 0; i_nei < k_nn; i_nei++){
            ind_nei = target_data[d_edg];
            position(i_nei+1,0) = float(weightk*(xyz_data[3 * ind_nei] - xyz_data[3 * i_ver])) ;
            position(i_nei+1,1) = float(weightk*(xyz_data[3 * ind_nei + 1] - xyz_data[3 * i_ver + 1])) ;
            position(i_nei+1,2) = float(weightk*(xyz_data[3 * ind_nei + 2] - xyz_data[3 * i_ver + 2])) ;
            d_edg++;
        }

        double four_norm = four_normal(xyz_data, target_data, i_ver, n_edg,nor_nei, xi_norm, weightk);

        float fn = (float)four_norm;
        
        position = fn * position;
        
        ei::MatrixXf centered_position = position.rowwise() - position.colwise().mean();
        ei::Matrix3f cov = (centered_position.adjoint() * centered_position) / float(k_nn + 1);
        ei::EigenSolver<Matrix3f> es(cov);
        //--- compute the eigen values and vectors---
        std::vector<float> ev = {es.eigenvalues()[0].real(),es.eigenvalues()[1].real(),es.eigenvalues()[2].real()};
        std::vector<int> indices(3);
        std::size_t n(0);
        std::generate(std::begin(indices), std::end(indices), [&]{ return n++; });
        std::sort(std::begin(indices),std::end(indices),
                       [&](int i1, int i2) { return ev[i1] > ev[i2]; } );
        std::vector<float> lambda = {(std::max(ev[indices[0]],0.f)),                        // 고유 벡터를 통해 정사영 했을 때, variance는 eigenvalue인데 var이 최대여야지 차원 감소시
                                    (std::max(ev[indices[1]],0.f)),                         // 원래의 형태를 잘 유지 할 수 있다.
                                    (std::max(ev[indices[2]],0.f))};
        ei::Vector3f v1 = {es.eigenvectors().col(indices[0])(0).real()                // 선형 변환의 주축
                               , es.eigenvectors().col(indices[0])(1).real()                // 선형 변환시 크기만 바뀌고 방향은 바뀌지 않는 벡터가 eigenvector
                               , es.eigenvectors().col(indices[0])(2).real()};              // 정사영 후 variance가 가장 큰 결과를 얻기 위해선 eigenvector에 정사영해야 한다.
        ei::Vector3f v2 = {es.eigenvectors().col(indices[1])(0).real()
                               , es.eigenvectors().col(indices[1])(1).real()
                               , es.eigenvectors().col(indices[1])(2).real()};
        ei::Vector3f v3 = {es.eigenvectors().col(indices[2])(0).real()
                               , es.eigenvectors().col(indices[2])(1).real()
                               , es.eigenvectors().col(indices[2])(2).real()};

        //--- compute the dimensionality features---
        float linearity  = (sqrtf(lambda[0]) - sqrtf(lambda[1])) / sqrtf(lambda[0]);
        float planarity  = (sqrtf(lambda[1]) - sqrtf(lambda[2])) / sqrtf(lambda[0]);
        float scattering =  sqrtf(lambda[2]) / sqrtf(lambda[0]);
        //--- compute the verticality feature---
        std::vector<float> unary_vector =
            {lambda[0] * fabsf(v1[0]) + lambda[1] * fabsf(v2[0]) + lambda[2] * fabsf(v3[0])
            ,lambda[0] * fabsf(v1[1]) + lambda[1] * fabsf(v2[1]) + lambda[2] * fabsf(v3[1])
            ,lambda[0] * fabsf(v1[2]) + lambda[1] * fabsf(v2[2]) + lambda[2] * fabsf(v3[2])};
        float norm = sqrt(unary_vector[0] * unary_vector[0] + unary_vector[1] * unary_vector[1]
                        + unary_vector[2] * unary_vector[2]);
        float verticality = unary_vector[2] / norm;

        float total_norm = sqrtf(linearity*linearity + planarity*planarity + scattering*scattering + verticality*verticality);
        
        ei::Vector3f vc1 = v2.cross(v3);
        ei::Vector3f vc2 = v3.cross(v1);
        ei::Vector3f vc3 = v1.cross(v2);
        
        ei::Matrix3f s;
        
        for(int i = 0 ; i <3 ; i++){
            s(0,i) = v1[i];
        }
        for(int q = 0 ; q <3 ; q++){
            s(1,q) = v2[q];
        }
        for(int k = 0 ; k <3 ; k++){
            s(2,k) = v3[k];
        }


        ei::EigenSolver<Matrix3f> esS(s);
        //--- compute the eigen values and vectors---
        std::vector<float> evS = {esS.eigenvalues()[0].real(),esS.eigenvalues()[1].real(),esS.eigenvalues()[2].real()};

        std::vector<int> tindiceS(3);
        std::size_t ns(0);
        std::generate(std::begin(tindiceS), std::end(tindiceS), [&]{ return ns++; });
        std::sort(std::begin(tindiceS),std::end(tindiceS),
                       [&](int j1, int j2) { return evS[j1] > evS[j2]; } );
        std::vector<float> lambdaS = {(std::max(evS[tindiceS[0]],0.f)),                       
                                    (std::max(evS[tindiceS[1]],0.f)),                     
                                    (std::max(evS[tindiceS[2]],0.f))};
        ei::Vector3f t1 = {esS.eigenvectors().col(tindiceS[0])(0).real()                
                            , esS.eigenvectors().col(tindiceS[0])(1).real()                
                            , esS.eigenvectors().col(tindiceS[0])(2).real()};            
        ei::Vector3f t2 = {esS.eigenvectors().col(tindiceS[1])(0).real()
                           , esS.eigenvectors().col(tindiceS[1])(1).real()
                           , esS.eigenvectors().col(tindiceS[1])(2).real()};
        ei::Vector3f t3 = {esS.eigenvectors().col(tindiceS[2])(0).real()
                           , esS.eigenvectors().col(tindiceS[2])(1).real()
                           , esS.eigenvectors().col(tindiceS[2])(2).real()};
        
        //---fill the geof vector---
        geof[i_ver][0] = linearity/total_norm;
        geof[i_ver][1] = planarity/total_norm;
        geof[i_ver][2] = scattering/total_norm;
        geof[i_ver][3] = verticality/total_norm;
        geof[i_ver][4] = float(t1[0]);
        geof[i_ver][5] = float(t1[1]);
        geof[i_ver][6] = float(t1[2]);

       
        //---progression---
        s_ver++;//if run in parellel s_ver behavior is udnefined, but gives a good indication of progress
        if (s_ver % 10000 == 0)
        {
            std::cout << s_ver << "% done          \r" << std::flush;
            std::cout << ceil(s_ver*100/n_ver) << "% done          \r" << std::flush;
        }
    }
    std::cout <<  std::endl;
    return VecvecToArray<float>::convert(geof);
}

PyObject * connected_comp(const uint32_t n_ver, const bpn::ndarray & source, const bpn::ndarray & target, const bpn::ndarray & active_edg, const int cutoff)
{//read data and run the L0-cut pursuit partition algorithm
    const uint32_t n_edg = bp::len(source);
    const uint32_t * source_data = reinterpret_cast<uint32_t*>(source.get_data());
    const uint32_t * target_data = reinterpret_cast<uint32_t*>(target.get_data());
    const char * active_edg_data = reinterpret_cast<char*>(active_edg.get_data());

    std::vector<uint32_t> in_component(n_ver,0);
    std::vector< std::vector<uint32_t> > components(1,std::vector<uint32_t>());

    connected_components(n_ver, n_edg, source_data, target_data, active_edg_data, in_component, components, cutoff);

    return to_py_tuple_components::convert(Components_tuple(components, in_component));
}


PyObject * random_subgraph(const int n_ver, const bpn::ndarray & source, const bpn::ndarray & target, const int subgraph_size)
{//read data and run the L0-cut pursuit partition algorithm

    const int n_edg = bp::len(source);
    const uint32_t * source_data = reinterpret_cast<uint32_t*>(source.get_data());
    const uint32_t * target_data = reinterpret_cast<uint32_t*>(target.get_data());

    std::vector<uint8_t> selected_edges(n_edg,0);
    std::vector<uint8_t> selected_vertices(n_ver,0);

    subgraph::random_subgraph(n_ver, n_edg, source_data, target_data, subgraph_size, selected_edges.data(), selected_vertices.data());

    return to_py_tuple_subgraph::convert(Subgraph_tuple(selected_edges,selected_vertices));
}

using namespace boost::python;
BOOST_PYTHON_MODULE(liblonet_c)
{
    _import_array();
    bp::to_python_converter<std::vector<std::vector<float>, std::allocator<std::vector<float> > >, VecvecToArray<float> >();
    bp::to_python_converter< Custom_tuple, to_py_tuple>();
    Py_Initialize();
    bpn::initialize();
    def("estimate_geof", estimate_geof);
    def("connected_comp", connected_comp);
    def("random_subgraph", random_subgraph);
}
