#ifndef OP_IMAGE_H
#define OP_IMAGE_H

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <string>
#include <vector>

#include "endian.h"
#include "utility.h"

 
namespace op {
    // Very simple class to write image files. Output only, and doesn't contain any image
    // processing tools (rotation, shapes, etc), only setting individual pixels. Image output is
    // not compressed.
    class Image {
    public:
        // Create a new image with given width and height.
        Image(int width, int height);

        // Set the pixel at (x, y) to the color (r, g, b, a).
        void set_pixel(int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a = 255);
        
        // Converts the image data to a PNG file and writes it as a sequence of char to *OutIter
        // incrementing it after every char.
        template<class OutIter>
        void write_png(OutIter& iter);
        
        // Returns the width of the image.
        const int width() const { return width_; }
        
        // Returns the height of the image.
        const int height() const { return height_; }

    private:
        int width_;
        int height_;

        std::vector<uint8_t> image_data;
    };
}



// Implementation.
namespace op {
    namespace detail {
        // Helper function to store a 32 bit big endian integer in memory.
        inline void store_be32(uint8_t* p, uint32_t x) {
            x = op::htobe(x);
            std::memcpy(p, reinterpret_cast<uint8_t*>(&x), sizeof(x));
        }


        // CRC implementation from PNG spec.
        inline uint32_t crc(uint8_t* buf, int len) {
            static std::array<uint32_t, 256> crc_table;
            static bool crc_table_computed = false;

            if (!crc_table_computed) {
                for (int n = 0; n < 256; ++n) {
                    uint32_t c = n;

                    for (int k = 0; k < 8; ++k) {
                        if (c & 1) c = 0xedb88320L ^ (c >> 1);
                        else c >>= 1;
                    }

                    crc_table[n] = c;
                }

                crc_table_computed = true;
            }
            
            uint32_t crc = 0xffffffff;
            for (int n = 0; n < len; ++n) {
                crc = crc_table[(crc ^ buf[n]) & 0xff] ^ (crc >> 8);
            }

            return crc ^ 0xffffffff;
        }


        // Adler32 checksum for deflate.
        inline uint32_t adler32(uint8_t* data, int len) {
            const int MOD_ADLER = 65521;
            uint32_t a = 1, b = 0;
         
            for (int i = 0; i < len; ++i) {
                a = (a + data[i]) % MOD_ADLER;
                b = (b + a) % MOD_ADLER;
            }
         
            return (b << 16) | a;
        }
    }


    inline Image::Image(int width, int height) {
        width_ = width;
        height_ = height;

        image_data.resize(width_ * height_ * 4);
    }


    inline void Image::set_pixel(int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a) {
        uint8_t *p = &image_data[4 * (width_ * y + x)];
        
        *p++ = r; *p++ = g; *p++ = b; *p++ = a;
    }


    template<class OutIter>
    inline void Image::write_png(OutIter& iter) {
        // PNG magic bytes.
        auto PNG_HEADER = make_array<uint8_t>(137, 80, 78, 71, 13, 10, 26, 10);
        iter.write(reinterpret_cast<const char*>(PNG_HEADER.data()), PNG_HEADER.size());

        // IHDR chunk.
        std::array<uint8_t, 4 + 4 + 13 + 4> ihdr_chunk;
        detail::store_be32(ihdr_chunk.data(), 13);                  // IHDR payload length.
        
        const std::string IHDR = "IHDR";
        std::copy(IHDR.begin(), IHDR.end(), ihdr_chunk.data() + 4); // IHDR identifier.

        detail::store_be32(ihdr_chunk.data() + 8,  width_);         // Width.
        detail::store_be32(ihdr_chunk.data() + 12, height_);        // Height.

        ihdr_chunk[16] = 8;                                         // Bitdepth.
        ihdr_chunk[17] = 6;                                         // Color type.
        ihdr_chunk[18] = 0;                                         // Compression type.
        ihdr_chunk[19] = 0;                                         // Filter type.
        ihdr_chunk[20] = 0;                                         // Interlace method.

        // Calculate CRC over everything but the payload length and write to stream.
        detail::store_be32(ihdr_chunk.data() + 21, detail::crc(ihdr_chunk.data() + 4, 17));
        iter.write(reinterpret_cast<const char*>(ihdr_chunk.data()), ihdr_chunk.size());

        // IDAT chunk.
        std::vector<uint8_t> idat_chunk(4); // reserve 4 bytes for the chunk size

        const std::string IDAT = "IDAT";
        std::copy(IDAT.begin(), IDAT.end(), std::back_inserter(idat_chunk)); // IDAT identifier

        // Every scanline must have a filter header byte.
        std::vector<uint8_t> scanlines;
        scanlines.reserve(4 * width_ * height_ + height_);

        for (int y = 0; y < height_; ++y) {
            scanlines.push_back(0); // no filter
            std::copy(image_data.begin() + 4 * width_ * y,
                      image_data.begin() + 4 * width_ * (y + 1),
                      std::back_inserter(scanlines));
        }

        // Deflate header.
        idat_chunk.push_back(120);
        idat_chunk.push_back(1);

        // Split in chunks of 65535.
        int num_chunks = 1 + scanlines.size() / 65535;
        for (int i = 0; i < scanlines.size(); ++i) {
            if (i % 65535 == 0) {
                --num_chunks;
                idat_chunk.push_back(num_chunks == 0); // last chunk?

                int chunk_size = 65535;
                if (num_chunks == 0) chunk_size = scanlines.size() % 65535;

                idat_chunk.push_back(chunk_size & 0xff);
                idat_chunk.push_back(chunk_size >> 8);
                idat_chunk.push_back(~chunk_size & 0xff);
                idat_chunk.push_back(~chunk_size >> 8);
            }

            idat_chunk.push_back(scanlines[i]);
        }

        // Add deflate checksum.
        idat_chunk.insert(idat_chunk.end(), 4, 0);
        detail::store_be32(idat_chunk.data() + idat_chunk.size() - 4,
                   detail::adler32(scanlines.data(), scanlines.size()));
        
        // Store payload size.
        int payload_size = idat_chunk.size() - 8;
        detail::store_be32(idat_chunk.data(), payload_size);

        // Store CRC and write.
        idat_chunk.insert(idat_chunk.end(), 4, 0);
        detail::store_be32(idat_chunk.data() + idat_chunk.size() - 4,
                   detail::crc(idat_chunk.data() + 4, payload_size + 4));
        iter.write(reinterpret_cast<const char*>(idat_chunk.data()), idat_chunk.size());

        // IEND chunk.
        auto IEND_CHUNK = make_array<uint8_t>(0, 0, 0, 0, 73, 69, 78, 68, 174, 66, 96, 130);
        iter.write(reinterpret_cast<const char*>(IEND_CHUNK.data()), IEND_CHUNK.size());
    }
}

#endif
