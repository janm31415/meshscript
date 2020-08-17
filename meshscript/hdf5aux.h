#pragma once

#include <hdf5.h>


template <class T>
struct get_hdf5_type
  {
  };

template <>
struct get_hdf5_type<float>
  {
  static hid_t get_type_id()
    {
    return H5T_IEEE_F32LE;
    }
  static hid_t get_mem_type_id()
    {
    return H5T_NATIVE_FLOAT;
    }
  };


template <>
struct get_hdf5_type<uint32_t>
  {
  static hid_t get_type_id()
    {
    return H5T_STD_U32LE;
    }
  static hid_t get_mem_type_id()
    {
    return H5T_NATIVE_UINT;
    }
  };

template <>
struct get_hdf5_type<uint16_t>
  {
  static hid_t get_type_id()
    {
    return H5T_STD_U16LE;
    }
  static hid_t get_mem_type_id()
    {
    return H5T_NATIVE_USHORT;
    }
  };

template <>
struct get_hdf5_type<int16_t>
  {
  static hid_t get_type_id()
    {
    return H5T_STD_I16LE;
    }
  static hid_t get_mem_type_id()
    {
    return H5T_NATIVE_SHORT;
    }
  };

template <>
struct get_hdf5_type<uint8_t>
  {
  static hid_t get_type_id()
    {
    return H5T_STD_U8LE;
    }
  static hid_t get_mem_type_id()
    {
    return H5T_NATIVE_UCHAR;
    }
  };


template <class T>
void save_attribute(hid_t group_id, const char* name, const T* data)
  {
  herr_t status;
  hsize_t dim = 1;
  auto dataspace_id = H5Screate_simple(1, &dim, NULL);
  auto attribute_id = H5Acreate2(group_id, name, get_hdf5_type<T>::get_type_id(), dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, get_hdf5_type<T>::get_mem_type_id(), data);
  status = H5Aclose(attribute_id);
  status = H5Sclose(dataspace_id);
  }

template <class T>
void save_2d_data(hid_t group_id, const char* name, size_t dim1, size_t dim2, const T* data)
  {
  herr_t status;
  hsize_t dims[2];
  dims[0] = dim1;
  dims[1] = dim2;
  auto dataspace_id = H5Screate_simple(2, dims, NULL);
  auto dataset_id = H5Dcreate2(group_id, name, get_hdf5_type<T>::get_type_id(), dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, get_hdf5_type<T>::get_mem_type_id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  }


template <class T>
void save_1d_data(hid_t group_id, const char* name, size_t dim1, const T* data)
  {
  herr_t status;
  hsize_t dim = dim1;
  auto dataspace_id = H5Screate_simple(1, &dim, NULL);
  auto dataset_id = H5Dcreate2(group_id, name, get_hdf5_type<T>::get_type_id(), dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, get_hdf5_type<T>::get_mem_type_id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  }


template <class T>
void load_attribute(hid_t group_id, const char* name, T* data)
  {
  herr_t status;
  hid_t attr = H5Aopen(group_id, name, H5P_DEFAULT);
  //auto space = H5Aget_space(attr);

  status = H5Aread(attr, get_hdf5_type<T>::get_mem_type_id(), data);

  status = H5Aclose(attr);
  //status = H5Sclose(space);  
  }


template <class T>
void load_1d_data(hid_t group_id, const char* name, size_t dim1, T* data)
  {
  herr_t status;
  auto dataset_id = H5Dopen2(group_id, name, H5P_DEFAULT);
  hsize_t dim = dim1;
  auto dataspace_id = H5Screate_simple(1, &dim, NULL);
  status = H5Dread(dataset_id, get_hdf5_type<T>::get_mem_type_id(), dataspace_id, H5S_ALL, H5P_DEFAULT, data);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  }


template <class T>
void load_2d_data(hid_t group_id, const char* name, size_t dim1, size_t dim2, T* data)
  {
  herr_t status;
  auto dataset_id = H5Dopen2(group_id, name, H5P_DEFAULT);
  hsize_t dims[2];
  dims[0] = dim1;
  dims[1] = dim2;
  auto dataspace_id = H5Screate_simple(2, dims, NULL);
  status = H5Dread(dataset_id, get_hdf5_type<T>::get_mem_type_id(), dataspace_id, H5S_ALL, H5P_DEFAULT, data);
  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  }