//
// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) Contributors to the OpenEXR Project.
//

#ifdef NDEBUG
#    undef NDEBUG
#endif

#include "testDeepScanLineBasic.h"
#include "random.h"

#include <assert.h>
#include <string.h>

#include <IlmThreadPool.h>
#include <ImfArray.h>
#include <ImfChannelList.h>
#include <ImfDeepFrameBuffer.h>
#include <ImfDeepScanLineInputFile.h>
#include <ImfDeepScanLineOutputFile.h>
#include <ImfPartType.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

namespace IMF = OPENEXR_IMF_NAMESPACE;
using namespace IMF;
using namespace std;
using namespace IMATH_NAMESPACE;
using namespace ILMTHREAD_NAMESPACE;

namespace
{

vector<int>           channelTypes;
Array2D<unsigned int> sampleCount;
Header                header;

void
generateRandomFile (
    const std::string filename,
    int               channelCount,
    Compression       compression,
    bool              bulkWrite,
    const Box2i&      dataWindow,
    const Box2i&      displayWindow)
{
    cout << "generating " << flush;
    header = Header (
        displayWindow,
        dataWindow,
        1,
        IMATH_NAMESPACE::V2f (0, 0),
        1,
        INCREASING_Y,
        compression);

    cout << "compression " << compression << " " << flush;

    int width  = dataWindow.max.x - dataWindow.min.x + 1;
    int height = dataWindow.max.y - dataWindow.min.y + 1;

    //
    // Add channels.
    //

    channelTypes.clear ();

    for (int i = 0; i < channelCount; i++)
    {
        int          type = random_int (3);
        stringstream ss;
        ss << i;
        string str = ss.str ();
        if (type == 0) header.channels ().insert (str, Channel (IMF::UINT));
        if (type == 1) header.channels ().insert (str, Channel (IMF::HALF));
        if (type == 2) header.channels ().insert (str, Channel (IMF::FLOAT));
        channelTypes.push_back (type);
    }

    header.setType (DEEPSCANLINE);


    //
    // bulkWrite mode stores the entire image;
    // otherwise only  single scanline is written
    //
    int bufferHeight = bulkWrite ? height : 1;

    Array<Array2D<void*>> data (channelCount);
    for (int i = 0; i < channelCount; i++)
        data[i].resizeErase (bufferHeight, width);

    sampleCount.resizeErase(height,width);

    //
    // for single scanline mode, a separate duplicate buffer with a single scanline
    //
    Array2D<unsigned int> localSampleCount(1,width);


    remove (filename.c_str ());
    DeepScanLineOutputFile file (filename.c_str (), header, 8);

    DeepFrameBuffer frameBuffer;

    if(bulkWrite)
    {
        frameBuffer.insertSampleCountSlice (Slice (
            IMF::UINT, // type // 7
            (char*) (&sampleCount[0][0] - dataWindow.min.x - dataWindow.min.y * width), // base
            sizeof (unsigned int) * 1,       // xStride
            sizeof (unsigned int) * width)); // yStride
    }
    else
    {
       frameBuffer.insertSampleCountSlice (Slice (
            IMF::UINT, // type // 7
            (char*) (&localSampleCount[0][0] - dataWindow.min.x), // base of single scanline
            sizeof (unsigned int) * 1,       // xStride
            0)); // yStride is zero, since all scanlines are stored in the same place
    }
    for (int i = 0; i < channelCount; i++)
    {
        PixelType type = NUM_PIXELTYPES;
        if (channelTypes[i] == 0) type = IMF::UINT;
        if (channelTypes[i] == 1) type = IMF::HALF;
        if (channelTypes[i] == 2) type = IMF::FLOAT;

        stringstream ss;
        ss << i;
        string str = ss.str ();

        int sampleSize = 0;
        if (channelTypes[i] == 0) sampleSize = sizeof (unsigned int);
        if (channelTypes[i] == 1) sampleSize = sizeof (half);
        if (channelTypes[i] == 2) sampleSize = sizeof (float);

        int pointerSize = sizeof (char*);

        if(bulkWrite)
        {
            frameBuffer.insert (
                str, // name // 6
                DeepSlice (
                    type, // type // 7
                    (char*) (&data[i][0][0] - dataWindow.min.x - dataWindow.min.y * width), // base // 8
                    pointerSize * 1,     // xStride// 9
                    pointerSize * width, // yStride// 10
                    sampleSize));        // sampleStride
        }
        else
        {
            frameBuffer.insert (
                str, // name // 6
                DeepSlice (
                    type, // type // 7
                    (char*) (&data[i][0][0] - dataWindow.min.x ), // base // 8
                    pointerSize * 1,     // xStride// 9
                    0, // yStride// 10
                    sampleSize));        // sampleStride
        }
    }

    file.setFrameBuffer (frameBuffer);

    int maxSamples = 10;

    bool bigFile = width > 1000 || height > 1000;

    cout << "writing " << flush;
    if (bulkWrite)
    {
        cout << "bulk " << flush;
        for (int i = 0; i < height; i++)
        {
            //
            // Fill in data at the last minute.
            //

            for (int j = 0; j < width; j++)
            {
                int samples = random_int (maxSamples);

                // big files write very sparse data for efficiency: most pixels have no samples

                if (bigFile && (i % 63 != 0 || j % 63 != 0)) { samples = 0; }
                sampleCount[i][j] = samples;

                for (int k = 0; k < channelCount; k++)
                {
                    if (samples > 0)
                    {

                        if (channelTypes[k] == 0)
                            data[k][i][j] = new unsigned int[sampleCount[i][j]];
                        if (channelTypes[k] == 1)
                            data[k][i][j] = new half[sampleCount[i][j]];
                        if (channelTypes[k] == 2)
                            data[k][i][j] = new float[sampleCount[i][j]];
                        for (unsigned int l = 0; l < sampleCount[i][j]; l++)
                        {
                            if (channelTypes[k] == 0)
                                ((unsigned int*) data[k][i][j])[l] =
                                    (i * width + j) % 2049;
                            if (channelTypes[k] == 1)
                                ((half*) data[k][i][j])[l] =
                                    (i * width + j) % 2049;
                            if (channelTypes[k] == 2)
                                ((float*) data[k][i][j])[l] =
                                    (i * width + j) % 2049;
                        }
                    }
                    else
                    {
                        data[k][i][j] = nullptr;
                    }
                }
            }
        }

        file.writePixels (height);

        for (int i = 0; i < height; i++)
          for (int j = 0; j < width; j++)
            for (int k = 0; k < channelCount; k++)
            {
                if (channelTypes[k] == 0) delete[](unsigned int*) data[k][i][j];
                if (channelTypes[k] == 1) delete[](half*) data[k][i][j];
                if (channelTypes[k] == 2) delete[](float*) data[k][i][j];
            }

    }
    else
    {
        cout << "per-line " << flush;
        for (int row = 0; row < height; row++)
        {
            //
            // Fill in data at the last minute.
            //

            for (int j = 0; j < width; j++)
            {
                int samples = random_int (maxSamples);

                // big files write very sparse data for efficiency: most pixels have no samples
                if (bigFile && (row % 63 != 0 || j % 63 != 0)) { samples = 0; }

                localSampleCount[0][j] = samples;
                sampleCount[row][j] = samples;

                for (int k = 0; k < channelCount; k++)
                {

                    if (samples > 0)
                    {
                        if (channelTypes[k] == 0)
                            data[k][0][j] = new unsigned int[localSampleCount[0][j]];
                        if (channelTypes[k] == 1)
                            data[k][0][j] = new half[localSampleCount[0][j]];
                        if (channelTypes[k] == 2)
                            data[k][0][j] = new float[localSampleCount[0][j]];
                        for (unsigned int l = 0; l < localSampleCount[0][j]; l++)
                        {
                            if (channelTypes[k] == 0)
                                ((unsigned int*) data[k][0][j])[l] =
                                    (row * width + j) % 2049;
                            if (channelTypes[k] == 1)
                                ((half*) data[k][0][j])[l] =
                                    (row * width + j) % 2049;
                            if (channelTypes[k] == 2)
                                ((float*) data[k][0][j])[l] =
                                    (row * width + j) % 2049;
                        }
                    }
                    else
                    {
                        data[k][0][j] = nullptr;
                    }
                }
            }
            file.writePixels (1);
            for (int j = 0; j < width; j++)
              for (int k = 0; k < channelCount; k++)
              {
                  if (channelTypes[k] == 0) delete[](unsigned int*) data[k][0][j];
                  if (channelTypes[k] == 1) delete[](half*) data[k][0][j];
                  if (channelTypes[k] == 2) delete[](float*) data[k][0][j];
              }
        }
    }


}

void
readFile (
    const std::string& filename,
    int                channelCount,
    bool               bulkRead,
    bool               randomChannels)
{
    if (randomChannels) { cout << " reading random channels " << flush; }
    else
    {
        cout << " reading all channels " << flush;
    }

    DeepScanLineInputFile file (filename.c_str (), 8);

    const Header& fileHeader = file.header ();
    assert (fileHeader.displayWindow () == header.displayWindow ());
    assert (fileHeader.dataWindow () == header.dataWindow ());
    assert (fileHeader.pixelAspectRatio () == header.pixelAspectRatio ());
    assert (fileHeader.screenWindowCenter () == header.screenWindowCenter ());
    assert (fileHeader.screenWindowWidth () == header.screenWindowWidth ());
    assert (fileHeader.lineOrder () == header.lineOrder ());
    assert (fileHeader.compression () == header.compression ());
    assert (fileHeader.channels () == header.channels ());
    assert (fileHeader.type () == header.type ());

    const Box2i& dataWindow = fileHeader.dataWindow ();

    int width  = dataWindow.max.x - dataWindow.min.x + 1;
    int height = dataWindow.max.y - dataWindow.min.y + 1;

    int bufferHeight = bulkRead ? height : 1;

    Array2D<unsigned int> localSampleCount;
    localSampleCount.resizeErase (bufferHeight, width);

    // also test filling channels. Generate up to 2 extra channels
    int fillChannels = random_int (3);

    Array<Array2D<void*>> data (channelCount + fillChannels);
    for (int i = 0; i < channelCount + fillChannels; i++)
        data[i].resizeErase (bufferHeight, width);

    DeepFrameBuffer frameBuffer;


    if(bulkRead)
    {
        frameBuffer.insertSampleCountSlice (Slice (
            IMF::UINT, // type // 7
            (char*) (&localSampleCount[0][0] - dataWindow.min.x - dataWindow.min.y * width), // base // 8)
            sizeof (unsigned int) * 1,       // xStride// 9
            sizeof (unsigned int) * width)); // yStride// 10
    }
    else
    {
            frameBuffer.insertSampleCountSlice (Slice (
            IMF::UINT, // type // 7
            (char*) (&localSampleCount[0][0] - dataWindow.min.x), // base // 8)
            sizeof (unsigned int) * 1,       // xStride// 9
            0)); // yStride// 10
    }

    vector<int> read_channel (channelCount);

    int channels_added = 0;

    for (int i = 0; i < channelCount; i++)
    {
        if (randomChannels) { read_channel[i] = random_int (2); }
        if (!randomChannels || read_channel[i] == 1)
        {
            PixelType type = NUM_PIXELTYPES;
            if (channelTypes[i] == 0) type = IMF::UINT;
            if (channelTypes[i] == 1) type = IMF::HALF;
            if (channelTypes[i] == 2) type = IMF::FLOAT;

            stringstream ss;
            ss << i;
            string str = ss.str ();

            int sampleSize = 0;
            if (channelTypes[i] == 0) sampleSize = sizeof (unsigned int);
            if (channelTypes[i] == 1) sampleSize = sizeof (half);
            if (channelTypes[i] == 2) sampleSize = sizeof (float);

            int pointerSize = sizeof (char*);

            if(bulkRead)
            {
                frameBuffer.insert (
                    str, // name // 6
                    DeepSlice (
                        type, // type // 7
                        (char*) (&data[i][0][0] - dataWindow.min.x - dataWindow.min.y * width), // base // 8)
                        pointerSize * 1,     // xStride// 9
                        pointerSize * width, // yStride// 10
                        sampleSize));        // sampleStride
            }
            else
            {
               frameBuffer.insert (
                    str, // name // 6
                    DeepSlice (
                        type, // type // 7
                        (char*) (&data[i][0][0] - dataWindow.min.x ), // base // 8)
                        pointerSize * 1,     // xStride// 9
                        0, // yStride// 10
                        sampleSize));        // sampleStride
            }
            channels_added++;
        }
    }

    if (channels_added == 0)
    {
        cout << "skipping " << flush;
        return;
    }
    for (int i = 0; i < fillChannels; ++i)
    {
        PixelType    type        = IMF::FLOAT;
        int          sampleSize  = sizeof (float);
        int          pointerSize = sizeof (char*);
        stringstream ss;
        // generate channel names that aren't in file but (might) interleave with existing file
        ss << i << "fill";
        string str = ss.str ();
        if(bulkRead)
        {
            frameBuffer.insert (str,                            // name // 6
                                    DeepSlice (type,                    // type // 7
                                    (char *) (&data[i+channelCount][0][0]
                                            - dataWindow.min.x
                                            - dataWindow.min.y * width),               // base // 8)
                                    pointerSize * 1,          // xStride// 9
                                    pointerSize * width,      // yStride// 10
                                    sampleSize));             // sampleStride
        }
        else
        {
            frameBuffer.insert (str,                            // name // 6
                                    DeepSlice (type,                    // type // 7
                                    (char *) (&data[i+channelCount][0][0]
                                            - dataWindow.min.x),               // base // 8)
                                    pointerSize * 1,          // xStride// 9
                                    0,      // yStride// 10
                                    sampleSize));             // sampleStride
        }
    }
    file.setFrameBuffer (frameBuffer);

    if (bulkRead)
    {
        cout << "bulk " << flush;
        file.readPixelSampleCounts (dataWindow.min.y, dataWindow.max.y);
        for (int i = 0; i < dataWindow.max.y - dataWindow.min.y + 1; i++)
        {
            for (int j = 0; j < width; j++)
                assert (localSampleCount[i][j] == sampleCount[i][j]);

            for (int j = 0; j < width; j++)
            {
                for (int k = 0; k < channelCount; k++)
                {
                    if (localSampleCount[i][j] > 0 &&
                        (!randomChannels || read_channel[k] == 1))
                    {
                        if (channelTypes[k] == 0)
                            data[k][i][j] =
                                new unsigned int[localSampleCount[i][j]];
                        if (channelTypes[k] == 1)
                            data[k][i][j] = new half[localSampleCount[i][j]];
                        if (channelTypes[k] == 2)
                            data[k][i][j] = new float[localSampleCount[i][j]];
                    }
                    else
                    {
                        data[k][i][j] = nullptr;
                    }
                }
                for (int f = 0; f < fillChannels; ++f)
                {
                    data[f + channelCount][i][j] =
                        new float[localSampleCount[i][j]];
                }
            }
        }

        file.readPixels (dataWindow.min.y, dataWindow.max.y);

        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
                for (int k = 0; k < channelCount; k++)
                {
                    if (!randomChannels || read_channel[k] == 1)
                    {
                        for (unsigned int l = 0; l < sampleCount[i][j]; l++)
                        {

                            if (channelTypes[k] == 0)
                            {
                                unsigned int* value =
                                    (unsigned int*) (data[k][i][j]);

                                if (value[l] !=
                                    static_cast<unsigned int> (i * width + j) %
                                        2049)
                                    cout << j << ", " << i << " error, should be "
                                        << (i * width + j) % 2049 << ", is "
                                        << value[l] << endl
                                        << flush;
                                assert (
                                    value[l] ==
                                    static_cast<unsigned int> (i * width + j) %
                                        2049);
                            }
                            if (channelTypes[k] == 1)
                            {
                                half* value = (half*) (data[k][i][j]);
                                if (value[l] != (i * width + j) % 2049)
                                    cout << j << ", " << i << " error, should be "
                                        << (i * width + j) % 2049 << ", is "
                                        << value[l] << endl
                                        << flush;
                                assert (
                                    ((half*) (data[k][i][j]))[l] ==
                                    (i * width + j) % 2049);
                            }
                            if (channelTypes[k] == 2)
                            {
                                float* value = (float*) (data[k][i][j]);
                                if (value[l] != (i * width + j) % 2049)
                                    cout << j << ", " << i << " error, should be "
                                        << (i * width + j) % 2049 << ", is "
                                        << value[l] << endl
                                        << flush;
                                assert (
                                    ((float*) (data[k][i][j]))[l] ==
                                    (i * width + j) % 2049);
                            }
                        }
                    }
                }

        for (int i = 0; i < height; i++)
            for (int j = 0; j < width; j++)
            {
                for (int k = 0; k < channelCount; k++)
                {
                    if (!randomChannels || read_channel[k] == 1)
                    {
                        if (channelTypes[k] == 0)
                            delete[](unsigned int*) data[k][i][j];
                        if (channelTypes[k] == 1) delete[](half*) data[k][i][j];
                        if (channelTypes[k] == 2) delete[](float*) data[k][i][j];
                    }
                }
                for (int f = 0; f < fillChannels; ++f)
                {
                    delete[](float*) data[f + channelCount][i][j];
                }
            }
    }

    else
    {
        cout << "per-line " << flush;
        for (int row = 0; row < dataWindow.max.y - dataWindow.min.y + 1; row++)
        {
            int y = row + dataWindow.min.y;
            file.readPixelSampleCounts (y);

            for (int j = 0; j < width; j++)
                assert (localSampleCount[0][j] == sampleCount[row][j]);

            for (int j = 0; j < width; j++)
            {
                for (int k = 0; k < channelCount; k++)
                {
                    if (localSampleCount[0][j] > 0 &&
                        (!randomChannels || read_channel[k] == 1))
                    {
                        if (channelTypes[k] == 0)
                            data[k][0][j] =
                                new unsigned int[localSampleCount[0][j]];
                        if (channelTypes[k] == 1)
                            data[k][0][j] = new half[localSampleCount[0][j]];
                        if (channelTypes[k] == 2)
                            data[k][0][j] = new float[localSampleCount[0][j]];
                    }
                    else
                    {
                        data[k][0][j] = nullptr;
                    }
                }
                for (int f = 0; f < fillChannels; ++f)
                {
                    data[f + channelCount][0][j] =
                        new float[localSampleCount[0][j]];
                }
            }

            file.readPixels (y);
            for (int j = 0; j < width; j++)
            {
                for (int k = 0; k < channelCount; k++)
                {
                    if (!randomChannels || read_channel[k] == 1)
                    {
                        for (unsigned int l = 0; l < sampleCount[row][j]; l++)
                        {
                            if (channelTypes[k] == 0)
                            {
                                unsigned int* value =
                                    (unsigned int*) (data[k][0][j]);
                                if (value[l] !=
                                    static_cast<unsigned int> (row * width + j) %
                                        2049)
                                    cout << j << ", " << row << " error, should be "
                                        << (row * width + j) % 2049 << ", is "
                                        << value[l] << endl
                                        << flush;
                                assert (
                                    value[l] ==
                                    static_cast<unsigned int> (row * width + j) %
                                        2049);
                            }
                            if (channelTypes[k] == 1)
                            {
                                half* value = (half*) (data[k][0][j]);
                                if (value[l] != (row * width + j) % 2049)
                                    cout << j << ", " << row << " error, should be "
                                        << (row * width + j) % 2049 << ", is "
                                        << value[l] << endl
                                        << flush;
                                assert (
                                    ((half*) (data[k][0][j]))[l] ==
                                    (row * width + j) % 2049);
                            }
                            if (channelTypes[k] == 2)
                            {
                                float* value = (float*) (data[k][0][j]);
                                if (value[l] != (row * width + j) % 2049)
                                    cout << j << ", " << row << " error, should be "
                                        << (row * width + j) % 2049 << ", is "
                                        << value[l] << endl
                                        << flush;
                                assert (
                                    ((float*) (data[k][0][j]))[l] ==
                                    (row * width + j) % 2049);
                            }
                        }
                    }
                }
            }

            for (int j = 0; j < width; j++)
            {
                for (int k = 0; k < channelCount; k++)
                {
                    if (!randomChannels || read_channel[k] == 1)
                    {
                        if (channelTypes[k] == 0)
                            delete[](unsigned int*) data[k][0][j];
                        if (channelTypes[k] == 1) delete[](half*) data[k][0][j];
                        if (channelTypes[k] == 2) delete[](float*) data[k][0][j];
                    }
                }
                for (int f = 0; f < fillChannels; ++f)
                {
                    delete[](float*) data[f + channelCount][0][j];
                }
            }
        }
    }
}

void
readWriteTest (
    const std::string& tempDir,
    int                channelCount,
    int                testTimes,
    const Box2i&       dataWindow,
    const Box2i&       displayWindow)
{
    cout << "Testing files with " << channelCount << " channels " << testTimes
         << " times." << endl
         << flush;

    std::string filename = tempDir + "imf_test_deep_scanline_basic.exr";

    for (int i = 0; i < testTimes; i++)
    {
        int         compressionIndex = i % 4;
        Compression compression;
        switch (compressionIndex)
        {
            case 0: compression = NO_COMPRESSION; break;
            case 1: compression = RLE_COMPRESSION; break;
            case 2: compression = ZIPS_COMPRESSION; break;
            case 3: compression = ZSTD_COMPRESSION; break;
        }

        generateRandomFile (
            filename,
            channelCount,
            compression,
            false,
            dataWindow,
            displayWindow);
        readFile (filename, channelCount, false, false);
        if (channelCount > 1) readFile (filename, channelCount, false, true);
        remove (filename.c_str ());
        cout << endl << flush;

        generateRandomFile (
            filename,
            channelCount,
            compression,
            true,
            dataWindow,
            displayWindow);
        readFile (filename, channelCount, true, false);
        if (channelCount > 1) readFile (filename, channelCount, true, true);
        remove (filename.c_str ());
        cout << endl << flush;
    }
}

void
testCompressionTypeChecks ()
{
    Header h;
    h.setType (DEEPTILE);
    h.channels ().insert ("Dummy", Channel ());
    h.compression () = NO_COMPRESSION;
    h.sanityCheck ();
    h.compression () = ZIPS_COMPRESSION;
    h.sanityCheck ();
    h.compression () = RLE_COMPRESSION;
    h.sanityCheck ();

    cout << "accepted valid compression types\n";
    //
    // these should fail
    //
    bool caught = false;
    try
    {
        h.compression () = ZIP_COMPRESSION;
        h.sanityCheck ();
        assert (false);
    }
    catch (...)
    {
        cout << "correctly identified bad compression setting (zip)\n";
        caught = true;
    }
    assert (caught);

    try
    {
        caught           = false;
        h.compression () = B44_COMPRESSION;
        h.sanityCheck ();
        assert (false);
    }
    catch (...)
    {
        cout << "correctly identified bad compression setting (b44)\n";
        caught = true;
    }
    assert (caught);

    try
    {
        caught           = false;
        h.compression () = B44A_COMPRESSION;
        h.sanityCheck ();
        assert (false);
    }
    catch (...)
    {
        cout << "correctly identified bad compression setting (b44a)\n";
        caught = true;
    }
    assert (caught);

    try
    {
        caught           = false;
        h.compression () = PXR24_COMPRESSION;
        h.sanityCheck ();
        assert (false);
    }
    catch (...)
    {
        cout << "correctly identified bad compression setting (pxr24)\n";
        caught = true;
    }
    assert (caught);
}

}; // namespace

namespace small
{
const int width  = 273;
const int height = 173;
const int minX   = 10;
const int minY   = 11;
} // namespace small

namespace large
{
const int width  = 5000;
const int height = 2500;
const int minX   = -22;
const int minY   = -12;
} // namespace large

void
testDeepScanLineBasic (const std::string& tempDir)
{
    try
    {
        cout << "\n\nTesting the DeepScanLineInput/OutputFile for basic use:\n"
             << endl;

        random_reseed (1);

        int numThreads = ThreadPool::globalThreadPool ().numThreads ();
        ThreadPool::globalThreadPool ().setNumThreads (4);

        testCompressionTypeChecks ();

        const Box2i largeDataWindow (
            V2i (large::minX, large::minY),
            V2i (
                large::minX + large::width - 1,
                large::minY + large::height - 1));
        const Box2i largeDisplayWindow (
            V2i (0, 0),
            V2i (
                large::minX + large::width * 2,
                large::minY + large::height * 2));

        readWriteTest (tempDir, 1, 3, largeDataWindow, largeDisplayWindow);

        const Box2i dataWindow (
            V2i (small::minX, small::minY),
            V2i (
                small::minX + small::width - 1,
                small::minY + small::height - 1));
        const Box2i displayWindow (
            V2i (0, 0),
            V2i (
                small::minX + small::width * 2,
                small::minY + small::height * 2));

        readWriteTest (tempDir, 1, 50, dataWindow, displayWindow);
        readWriteTest (tempDir, 3, 25, dataWindow, displayWindow);
        readWriteTest (tempDir, 10, 10, dataWindow, displayWindow);

        ThreadPool::globalThreadPool ().setNumThreads (numThreads);

        cout << "ok\n" << endl;
    }
    catch (const std::exception& e)
    {
        cerr << "ERROR -- caught exception: " << e.what () << endl;
        assert (false);
    }
}
