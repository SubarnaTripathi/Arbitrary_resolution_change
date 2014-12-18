Arbitrary_resolution_change
===========================

This can change the resolution of any 420 yuv image/video to any arbitrary user input resolution. Among the three input algorithms supported (nearest neighbor, bilinear filter and polyphase filter), the polyphase filter implementation reuses a significant code from ffmpeg library. I only simplified some parts.

Example use:

user would be asked to enter yuv file name, then height, width of the original image (should be multiple of 16)
User will also be asked for output resolution and filter type and the frame number.
filter type = 1 => nearest neighbor
filter type = 2 => bilinear filter
filter type = 3 => polyphase filter
