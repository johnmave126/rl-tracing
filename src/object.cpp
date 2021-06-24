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

#include <tracer/object.h>

TRACER_NAMESPACE_BEGIN

void TracerObject::addChild(TracerObject *) {
    throw TracerException(
        "TracerObject::addChild() is not implemented for objects of type '%s'!",
        classTypeName(getClassType()));
}

void TracerObject::activate() { /* Do nothing */ }
void TracerObject::setParent(TracerObject *) { /* Do nothing */ }

std::map<std::string, TracerObjectFactory::Constructor> *TracerObjectFactory::m_constructors = nullptr;

void TracerObjectFactory::registerClass(const std::string &name, const Constructor &constr) {
    if (!m_constructors)
        m_constructors = new std::map<std::string, TracerObjectFactory::Constructor>();
    (*m_constructors)[name] = constr;
}

TRACER_NAMESPACE_END
