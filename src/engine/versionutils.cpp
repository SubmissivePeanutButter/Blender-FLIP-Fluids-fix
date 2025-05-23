/*
MIT License

Copyright (C) 2025 Ryan L. Guy & Dennis Fassbaender

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "versionutils.h"

int VersionUtils::_major = 1;
int VersionUtils::_minor = 8;
int VersionUtils::_revision = 3;
std::string VersionUtils::_label = "1.8.3 GitHub Release 23-MAR-2025";

void VersionUtils::getVersion(int *major, int *minor, int *revision) {
    *major = _major;
    *minor = _minor;
    *revision = _revision;
}

int VersionUtils::getMajor() {
    return _major;
}

int VersionUtils::getMinor() {
    return _minor;
}

int VersionUtils::getRevision() {
    return _revision;
}

std::string VersionUtils::getLabel() {
    return _label;
}
