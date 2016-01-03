#ifndef image_h
#define image_h

/* Image type - contains height, width, and data */
struct texImage {
    unsigned long sizeX;
    unsigned long sizeY;
    char *data;
};

int ImageLoad(char *filename, texImage *image);


#endif
