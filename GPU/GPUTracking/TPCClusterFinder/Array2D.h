// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PackedCharge.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_ARRAY2D_H
#define O2_GPU_ARRAY2D_H

#include "clusterFinderDefs.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

template <typename T, size_t Width, size_t Height>
class TilingLayoutArray2D
{
 public:
  GPUdi() explicit TilingLayoutArray2D(GPUglobalref() T* d) : data(d) {}

  GPUdi() T& operator[](const ChargePos& p) { return data[idx(p)]; }
  GPUdi() const T& operator[](const ChargePos& p) const { return data[idx(p)]; }

 private:
  GPUglobalref() T* data;

  GPUdi() size_t idx(const ChargePos& p) const
  {
    size_t time = p.time + PADDING_TIME;

    const size_t widthInTiles = (TPC_NUM_OF_PADS + Width - 1) / Width;

    const size_t tilePad = p.gpad / Width;
    const size_t tileTime = time / Height;

    const size_t inTilePad = p.gpad % Width;
    const size_t inTileTime = time % Height;

    return (tileTime * widthInTiles + tilePad) * (Width * Height) + inTileTime * Width + inTilePad;
  }
};

template <typename T>
class LinearLayoutArray2D
{
 public:
  GPUdi() explicit LinearLayoutArray2D(GPUglobalref() T* d) : data(d) {}

  GPUdi() T& operator[](const ChargePos& p) { return data[idx(p)]; }
  GPUdi() const T& operator[](const ChargePos& p) const { return data[idx(p)]; }

 private:
  GPUglobalref() T* data;

  GPUdi() size_t idx(const ChargePos& p) const
  {
    size_t time = p.time + PADDING_TIME;
    return TPC_NUM_OF_PADS * time + p.gpad;
  }
};

template <typename T>
class Array2DMapper;

#if defined(CHARGEMAP_TILING_LAYOUT)
template <>
class Array2DMapper<PackedCharge>
{
 public:
  using Type = TilingLayoutArray2D<PackedCharge, 4, 4>;
};

template <>
class Array2DMapper<uchar>
{
 public:
  using Type = TilingLayoutArray2D<uchar, 8, 8>;
};
#else
template <typename T>
class Array2DMapper
{
 public:
  using Type = LinearLayoutArray2D<T>;
};
#endif

template <typename T>
using Array2D = typename Array2DMapper<T>::Type;

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#define CHARGE(map, gpad, time) map[{gpad, time}]
#define IS_PEAK(map, gpad, time) map[{gpad, time}]

#endif
