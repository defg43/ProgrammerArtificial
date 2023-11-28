#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "genann.h"

#define lengthof(x) (sizeof(x)/sizeof(x[0]))

typedef struct {
    union {
        unsigned char channel[4];
        struct {
            unsigned char r, g, b, a;
        };
    };
} pixel_t;

typedef struct {
    union {
        double channel[4];
        struct {
            double r, g, b, a;
        };
    };
} double_pixel_t;

typedef enum : uint64_t {
    wood = 1,
    stone,
    ores,
    stones,
    workbenches
} texture_tag_t;

typedef struct {
    texture_tag_t tag;
    double_pixel_t data[16][16];
} texture_t;

// Function to load a 16x16 texture from a file
int loadTextureFromFile(const char* filename, texture_t* texture) {
    int width, height, channels;
    unsigned char* image = stbi_load(filename, &width, &height, &channels, 4); // 4 channels for RGBA

    if (!image) {
        // Error loading the image
        return 0;
    }

    if (width != 16 || height != 16) {
        // Ensure the image is 16x16
        stbi_image_free(image);
        return 0;
    }

    texture->tag = wood; // Set the appropriate texture tag

    // Copy pixel data
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            int index = (i * 16 + j) * 4; // Each pixel has 4 channels (RGBA)
            texture->data[i][j].r = image[index] / 255.0;
            texture->data[i][j].g = image[index + 1] / 255.0;
            texture->data[i][j].b = image[index + 2] / 255.0;
            texture->data[i][j].a = image[index + 3] / 255.0;
        }
    }

    stbi_image_free(image);
    return 1;
}


void saveTextureToFile(const char *filename, const texture_t *generatedTexture) {
    unsigned char *imageData = (unsigned char *)malloc(4 * 16 * 16);
    if (!imageData) {
        fprintf(stderr, "Error: Unable to allocate memory for image data.\n");
        return;
    }
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            int index = (i * 16 + j) * 4;
            imageData[index] = (unsigned char)(generatedTexture->data[i][j].r * 255.0);
            imageData[index + 1] = (unsigned char)(generatedTexture->data[i][j].g * 255.0);
            imageData[index + 2] = (unsigned char)(generatedTexture->data[i][j].b * 255.0);
            imageData[index + 3] = (unsigned char)(generatedTexture->data[i][j].a * 255.0);
        }
    }
    if (stbi_write_png(filename, 16, 16, 4, imageData, 4 * 16) == 0) {
        fprintf(stderr, "Error: Unable to save the image to %s\n", filename);
    } else {
        printf("Generated image saved to %s\n", filename);
    }
    free(imageData);
}

void genann_visualize(genann const *ann, const char *filename) {
    // Create a 2D array to represent the image (assuming square image for simplicity).
    const int image_size = 256;  // Adjust as needed.
    unsigned char image[image_size][image_size][3];  // RGB format.

    // Initialize the image to white.
    memset(image, 255, sizeof(image));

    // Calculate the position of neurons in the image.
    int input_offset = image_size / (ann->inputs + 1);
    int hidden_offset = image_size / (ann->hidden * ann->hidden_layers + 1);
    int output_offset = image_size / (ann->outputs + 1);

    // Draw neurons on the image.
    for (int i = 0; i < ann->inputs; ++i) {
        image[input_offset * (i + 1)][50][0] = 0;  // Red
        image[input_offset * (i + 1)][50][1] = 0;  // Green
        image[input_offset * (i + 1)][50][2] = 0;  // Blue
    }

    for (int i = 0; i < ann->hidden_layers; ++i) {
        for (int j = 0; j < ann->hidden; ++j) {
            image[hidden_offset * (i + 1)][j * hidden_offset][0] = 0;
            image[hidden_offset * (i + 1)][j * hidden_offset][1] = 0;
            image[hidden_offset * (i + 1)][j * hidden_offset][2] = 0;
        }
    }

    for (int i = 0; i < ann->outputs; ++i) {
        image[image_size - output_offset * (i + 1)][50][0] = 0;
        image[image_size - output_offset * (i + 1)][50][1] = 0;
        image[image_size - output_offset * (i + 1)][50][2] = 0;
    }

    // Draw connections between neurons.
    for (int i = 0; i < ann->inputs; ++i) {
        for (int j = 0; j < ann->hidden; ++j) {
            int x1 = input_offset * (i + 1);
            int y1 = 50;
            int x2 = hidden_offset;
            int y2 = j * hidden_offset;
            for (int k = 0; k < image_size / 2; ++k) {
                image[x1 + k * (x2 - x1) / (image_size / 2)][y1 + k * (y2 - y1) / (image_size / 2)][0] = 0;
                image[x1 + k * (x2 - x1) / (image_size / 2)][y1 + k * (y2 - y1) / (image_size / 2)][1] = 0;
                image[x1 + k * (x2 - x1) / (image_size / 2)][y1 + k * (y2 - y1) / (image_size / 2)][2] = 0;
            }
        }
    }

    for (int i = 0; i < ann->hidden_layers - 1; ++i) {
        for (int j = 0; j < ann->hidden; ++j) {
            for (int k = 0; k < ann->hidden; ++k) {
                int x1 = hidden_offset * (i + 1);
                int y1 = j * hidden_offset;
                int x2 = hidden_offset * (i + 2);
                int y2 = k * hidden_offset;
                for (int l = 0; l < image_size / 2; ++l) {
                    image[x1 + l * (x2 - x1) / (image_size / 2)][y1 + l * (y2 - y1) / (image_size / 2)][0] = 0;
                    image[x1 + l * (x2 - x1) / (image_size / 2)][y1 + l * (y2 - y1) / (image_size / 2)][1] = 0;
                    image[x1 + l * (x2 - x1) / (image_size / 2)][y1 + l * (y2 - y1) / (image_size / 2)][2] = 0;
                }
            }
        }
    }

    for (int i = 0; i < ann->hidden; ++i) {
        for (int j = 0; j < ann->outputs; ++j) {
            int x1 = hidden_offset * ann->hidden_layers;
            int y1 = i * hidden_offset;
            int x2 = image_size - output_offset * (j + 1);
            int y2 = 50;
            for (int k = 0; k < image_size / 2; ++k) {
                image[x1 + k * (x2 - x1) / (image_size / 2)][y1 + k * (y2 - y1) / (image_size / 2)][0] = 0;
                image[x1 + k * (x2 - x1) / (image_size / 2)][y1 + k * (y2 - y1) / (image_size / 2)][1] = 0;
                image[x1 + k * (x2 - x1) / (image_size / 2)][y1 + k * (y2 - y1) / (image_size / 2)][2] = 0;
            }
        }
    }

    // Save the image to a PPM file.
    FILE *fp = fopen(filename, "wb");
    fprintf(fp, "P6\n%d %d\n255\n", image_size, image_size);
    fwrite(image, sizeof(image), 1, fp);
    fclose(fp);
}


void visualizeWeightsAsImage(genann const *ann, const char *filename) {
    const int width = ann->total_weights / 3;  // Assuming RGB channels for simplicity.
    const int height = 1;

    unsigned char *image = (unsigned char *)malloc(width * height * 3);  // RGB format

    if (!image) {
        fprintf(stderr, "Error: Unable to allocate memory for image data.\n");
        return;
    }

    // Normalize the weights to the range [0, 255].
    double minWeight = ann->weight[0];
    double maxWeight = ann->weight[0];

    for (int i = 0; i < ann->total_weights; ++i) {
        if (ann->weight[i] < minWeight) {
            minWeight = ann->weight[i];
        }
        if (ann->weight[i] > maxWeight) {
            maxWeight = ann->weight[i];
        }
    }

    // Map weights to RGB values and set pixels in the image.
    for (int i = 0; i < width; ++i) {
        double normalizedWeight = (ann->weight[i] - minWeight) / (maxWeight - minWeight);
        unsigned char pixelValue = (unsigned char)(normalizedWeight * 255.0);
        image[i * 3] = pixelValue;      // Red channel
        image[i * 3 + 1] = pixelValue;  // Green channel
        image[i * 3 + 2] = pixelValue;  // Blue channel
    }

    // Save the image to a file.
    if (stbi_write_png(filename, width, height, 3, image, width * 3) == 0) {
        fprintf(stderr, "Error: Unable to save the image to %s\n", filename);
    } else {
        printf("Weight visualization image saved to %s\n", filename);
    }

    free(image);
}

void generateImageFromAnn(genann const *ann, double const *input, const char *filename) {
    const int imageSize = 16;  // Assuming a 16x16 image.
    unsigned char image[imageSize][imageSize][4];  // 4 channels (RGBA) format.

    // Run the neural network on the input.
    double const *output = genann_run(ann, input);

    // Map the output to RGBA values and set pixels in the image.
    for (int i = 0; i < imageSize; ++i) {
        for (int j = 0; j < imageSize; ++j) {
            int index = (i * imageSize + j) * 4;  // Each pixel has 4 channels (RGBA).
            image[i][j][0] = (unsigned char)(output[index] * 255.0);    // Red channel
            image[i][j][1] = (unsigned char)(output[index + 1] * 255.0);  // Green channel
            image[i][j][2] = (unsigned char)(output[index + 2] * 255.0);  // Blue channel
            image[i][j][3] = (unsigned char)(output[index + 3] * 255.0);  // Alpha channel
        }
    }

    // Save the image to a file.
    if (stbi_write_png(filename, imageSize, imageSize, 4, image, imageSize * 4) == 0) {
        fprintf(stderr, "Error: Unable to save the image to %s\n", filename);
    } else {
        printf("Image generated from neural network saved to %s\n", filename);
    }
}

int main() {
    char *filename[] = {
        "./oak_planks.new.png",
        "./oak_planks.old.png",
    };
    texture_t texture[lengthof(filename)];
    printf("%ld\n", lengthof(filename));
    for (size_t i = 0; i < lengthof(filename); i++) {
        if (loadTextureFromFile(filename[i], &texture[i])) {
            printf("successfully loaded %s\n", filename[i]);
        } else {
            printf("failed to load %s\n", filename[i]);
        }
    }

    double input[4 * 16][16];
    memcpy(input, texture[0].data, sizeof(input));

    double output[4 * 16][16];
    memcpy(output, texture[1].data, sizeof(output));

    genann *ann = genann_init(16 * 16 * 4, 16, 16 * 4, 16 * 16 * 4);

	printf("we do a bit of training ...\n");
    for (int iterations = 0; iterations < 100; iterations++) {
        	for (unsigned int i = 0; i < lengthof(input); i++)
            	genann_train(ann, input[i], output[i], 3);
    }
	printf("training should be done\n");
	genann_visualize(ann, "visualisation.png");
	generateImageFromAnn(ann, *input, "result.png");
	genann_free(ann);
    return 0;
}
