/*
    This file is part of Tracer, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    [redacted] is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    [redacted] is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <tracer/color.h>
#include <tracer/vector.h>

TRACER_NAMESPACE_BEGIN

/**
 * \brief Stores a RGB high dynamic-range bitmap
 *
 * The bitmap class provides I/O support using the OpenEXR file format
 */
class Bitmap : public Eigen::Array<Color3f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> {
public:
    typedef Eigen::Array<Color3f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Base;

    /**
     * \brief Allocate a new bitmap of the specified size
     *
     * The contents will initially be undefined, so make sure
     * to call \ref clear() if necessary
     */
    Bitmap(const Vector2i &size = Vector2i(0, 0))
        : Base(size.y(), size.x()) { }

    /// Load an OpenEXR file with the specified filename
    Bitmap(const std::string &filename);

    /// Save the bitmap as an EXR file with the specified filename
    void save(const std::string &filename);
};

TRACER_NAMESPACE_END
