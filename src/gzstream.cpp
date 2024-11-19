// ============================================================================
// gzstream, C++ iostream classes wrapping the zlib compression library.
// Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : gzstream.C
// Revision      : $Revision: 1.7 $
// Revision_date : $Date: 2003/01/08 14:41:27 $
// Author(s)     : Deepak Bandyopadhyay, Lutz Kettner
// 
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library".
// ============================================================================

#include "gzstream.h"
#include <iostream>
#include <string.h>  // for memcpy
#include <stdexcept>

#include <iostream>

#ifdef GZSTREAM_NAMESPACE
namespace GZSTREAM_NAMESPACE {
#endif

// ----------------------------------------------------------------------------
// Internal classes to implement gzstream. See header file for user classes.
// ----------------------------------------------------------------------------

// --------------------------------------
// class gzstreambuf:
// --------------------------------------

void get_fmode( char* fmode, int open_mode) {
    char* fmodeptr = fmode;
    if (open_mode & std::ios::in)
        *fmodeptr++ = 'r';
    *fmodeptr++ = 'b';
    *fmodeptr = '\0';
}

gzstreambuf* gzstreambuf::open( const char* name, int open_mode) {
    filename = std::string(name);
    // const size_t n = strlen(name);  // excludes null terminator
    // filename = new char[n + 1]{}; // {} zero initializes array
    // std::copy_n(name, n, filename);
    
    if ( is_open())
        return (gzstreambuf*)0;
    mode = open_mode;
    // no append nor read/write mode
    if ((mode & std::ios::ate) || (mode & std::ios::app)
        || ((mode & std::ios::in) && (mode & std::ios::out)))
        return (gzstreambuf*)0;
    get_fmode(fmode, mode);
    handle = fopen(filename.c_str(), fmode);
    if (!handle) {
        return 0;
    }
    
    file = gzdopen(fileno(handle), fmode);
    if (file == 0)
        return (gzstreambuf*)0;
    opened = 1;
    return this;
}

gzstreambuf * gzstreambuf::close() {
    if ( is_open()) {
        sync();
        opened = 0;
        if ( gzclose( file) == Z_OK)
            return this;
    }
    return (gzstreambuf*)0;
}

/// @brief seek the bgzip file to a new offset
/// @param offset file offset of new bgzip block to seek to
void gzstreambuf::seek(bcf::Offsets offset) {
    // reset the input buffer
    setp( buffer, buffer + (bufferSize-1));
    setg( buffer + 4,  // beginning of putback area
          buffer + 4,  // read position
          buffer + 4); // end position
    
    // close old gzfile object and associated file handle
    fclose(handle);
    gzclose(file);

    // open new file handle, and seek to an offset in the compressed file
    handle = fopen(filename.c_str(), fmode);
    fseek(handle, offset.c_offset, SEEK_SET);
    ftell(handle);
    
    // open new gzfile object using the file descriptor (at new offset);
    file = gzdopen(fileno(handle), fmode);

    if (file == 0) {
        throw std::invalid_argument("cannot seek within this gzfile");
    }
    
    // read through to the correct offset in the uncompressed data. The 
    // uncompressed offset must be less than 2 ** 16, since that is the max
    // BGZF chunk size (uncompressed or compressed).
    char tmp[65536];
    gzread(file, tmp, offset.u_offset);
}

int gzstreambuf::underflow() { // used for input buffer only
    if ( gptr() && ( gptr() < egptr()))
        return * reinterpret_cast<unsigned char *>( gptr());

    if ( ! (mode & std::ios::in) || ! opened)
        return EOF;
    // Josuttis' implementation of inbuf
    int n_putback = gptr() - eback();
    if ( n_putback > 4)
        n_putback = 4;
    memcpy( buffer + (4 - n_putback), gptr() - n_putback, n_putback);

    int num = gzread( file, buffer+4, bufferSize-4);
    if (num <= 0) // ERROR or EOF
        return EOF;

    // reset buffer pointers
    setg( buffer + (4 - n_putback),   // beginning of putback area
          buffer + 4,                 // read position
          buffer + 4 + num);          // end of buffer

    // return next character
    return * reinterpret_cast<unsigned char *>( gptr());    
}

// --------------------------------------
// class gzstreambase:
// --------------------------------------

gzstreambase::gzstreambase( const char* name, int mode) {
    init( &buf);
    open( name, mode);
}

gzstreambase::~gzstreambase() {
    buf.close();
}

void gzstreambase::open( const char* name, int open_mode) {
    if ( ! buf.open( name, open_mode)) {
        clear( rdstate() | std::ios::badbit);
    }
}

void gzstreambase::close() {
    if ( buf.is_open()) {
        if ( ! buf.close()) {
            clear( rdstate() | std::ios::badbit);
        }
    }
}

#ifdef GZSTREAM_NAMESPACE
} // namespace GZSTREAM_NAMESPACE
#endif

// ============================================================================
// EOF //
