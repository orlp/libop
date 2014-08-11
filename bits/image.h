#ifndef OP_IMAGE_H
#define OP_IMAGE_H

#include <iostream>
#include <cstdint>
#include <array>

#include "endian.h"

 
namespace op {
    class Image {
    public:
        Image(int width, int height) {
            width_ = width;
            height_ = height;

            image_data.resize(width_ * height_ * 4);
        }

        void set_pixel(int x, int y, unsigned char red, unsigned char green, unsigned char blue, unsigned char alpha = 255) {
            unsigned char *p = &image_data[4 * (width_ * y + x)];
            
            *p++ = red;
            *p++ = green;
            *p++ = blue;
            *p++ = alpha;
        }

        template<class OutIter>
        void write_png(OutIter& iter) {
            // PNG magic bytes
            const std::array<unsigned char, 8> PNG_HEADER = { 137, 80, 78, 71, 13, 10, 26, 10 };
            iter.write(reinterpret_cast<const char*>(PNG_HEADER.data()), PNG_HEADER.size());

            // IHDR chunk
            std::array<unsigned char, 4 + 4 + 13 + 4> ihdr_chunk;
            store_be32(ihdr_chunk.data(), 13);                          // IHDR payload length
            
            const std::string IHDR = "IHDR";
            std::copy(IHDR.begin(), IHDR.end(), ihdr_chunk.data() + 4); // IHDR identifier

            store_be32(ihdr_chunk.data() + 8,  width_);                 // width
            store_be32(ihdr_chunk.data() + 12, height_);                // height

            ihdr_chunk[16] = 8;                                         // bitdepth
            ihdr_chunk[17] = 6;                                         // color type (truecolor with alpha)
            ihdr_chunk[18] = 0;                                         // compression type
            ihdr_chunk[19] = 0;                                         // filter type
            ihdr_chunk[20] = 0;                                         // interlace method

            // calculate CRC over everything but the payload length and write to stream
            store_be32(ihdr_chunk.data() + 21, crc(ihdr_chunk.data() + 4, 17));
            iter.write(reinterpret_cast<const char*>(ihdr_chunk.data()), ihdr_chunk.size());

            // IDAT chunk
            std::vector<unsigned char> idat_chunk(4); // reserve 4 bytes for the chunk size

            const std::string IDAT = "IDAT";
            std::copy(IDAT.begin(), IDAT.end(), std::back_inserter(idat_chunk)); // IDAT identifier

            // every scanline must have a filter header byte
            std::vector<unsigned char> scanlines;
            scanlines.reserve(4 * width_ * height_ + height_);

            for (int y = 0; y < height_; ++y) {
                scanlines.push_back(0); // no filter
                std::copy(image_data.begin() + 4 * width_ * y, image_data.begin() + 4 * width_ * (y + 1), std::back_inserter(scanlines));
            }

            // deflate header
            idat_chunk.push_back(120);
            idat_chunk.push_back(1);

            // split in chunks of 65535
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

            // add deflate checksum
            idat_chunk.insert(idat_chunk.end(), 4, 0);
            store_be32(idat_chunk.data() + idat_chunk.size() - 4, adler32(scanlines.data(), scanlines.size()));
            
            // store payload size
            int payload_size = idat_chunk.size() - 8;
            store_be32(idat_chunk.data(), payload_size);

            // store CRC and write
            idat_chunk.insert(idat_chunk.end(), 4, 0);
            store_be32(idat_chunk.data() + idat_chunk.size() - 4, crc(idat_chunk.data() + 4, payload_size + 4));
            iter.write(reinterpret_cast<const char*>(idat_chunk.data()), idat_chunk.size());

            // IEND chunk
            const std::array<unsigned char, 12> IEND_CHUNK = {0, 0, 0, 0, 73, 69, 78, 68, 174, 66, 96, 130};
            iter.write(reinterpret_cast<const char*>(IEND_CHUNK.data()), IEND_CHUNK.size());
        }

    private:
        int width_;
        int height_;

        std::vector<unsigned char> image_data;

        void store_be32(unsigned char* p, uint32_t x) {
            x = op::endian::htobe(x);
            std::memcpy(p, reinterpret_cast<unsigned char*>(&x), sizeof(x));
        }

        // CRC implementation from PNG spec
        uint32_t crc(unsigned char *buf, int len) {
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

        // adler32 checksum for deflate
        uint32_t adler32(unsigned char *data, int32_t len) {
            const int MOD_ADLER = 65521;
            uint32_t a = 1, b = 0;
         
            for (int i = 0; i < len; ++i) {
                a = (a + data[i]) % MOD_ADLER;
                b = (b + a) % MOD_ADLER;
            }
         
            return (b << 16) | a;
        }
    };
}

#endif
