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

#include <tracer/common.h>
#include <nanogui/screen.h>

TRACER_NAMESPACE_BEGIN

class TracerScreen : public nanogui::Screen {
public:
    TracerScreen(const ImageBlock &block);
    virtual ~TracerScreen();

    void drawContents();
private:
    const ImageBlock &m_block;
    nanogui::GLShader *m_shader = nullptr;
    nanogui::Slider *m_slider = nullptr;
    uint32_t m_texture = 0;
    float m_scale = 1.f;
};

TRACER_NAMESPACE_END
