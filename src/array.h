namespace array {
    float* rand_array(size_t n_elements)
    {
        float* array = (float *)malloc(n_elements * sizeof(float));
        for (size_t i = 0; i < n_elements; i++) {
            array[i] = udist(rng);
        }

        return array;
    }

    void fill_array(float *array, float value, size_t n_elements)
    {
        for(int i = 0; i < n_elements; i++){
            array[i] = value;
        }
    }


    struct ResizableArray
    {
        size_t max_points;
        size_t cur_points;
        float *data;
    };

    ResizableArray create_resizable_array(size_t max_points) {
        ResizableArray array{max_points, 0, (float *)malloc(3 * max_points * sizeof(float))};
        return array;
    }

    void resize_array(ResizableArray *array, size_t new_points)
    {
        array->max_points = new_points;
        array->data = (float *)realloc(array->data, 3 * new_points * sizeof(float));
    }

    void fill_resizable_array(ResizableArray *array, size_t start, size_t end,
                            float fill_value)
    {
        for (; start < end; start++)
        {
            array->data[start] = fill_value;
        }
    }

    ResizableArray* create_empty_data_structure(int n_processors, size_t max_points) {
        ResizableArray* data_structure = (ResizableArray*) malloc(n_processors*sizeof(ResizableArray));
        for (size_t i = 0; i < n_processors; i++) {
            data_structure[i] = create_resizable_array(max_points);
        }
    }

    void insert_point(ResizableArray *DS, float x, float y, float z, float R, int n_processors) {
        float step = 1.f / (float)n_processors;
        float low = 0.f;
        float high = low + step;

        for(size_t j = 0; j < n_processors; j++) {
            if (x >= low-R && x < high+R) {
                if (DS[j].max_points < DS[j].cur_points + 1) {
                    resize_array(&DS[j], DS[j].max_points * 2);
                }

                DS[j].cur_points += 1;
                DS[j].data[3*DS[j].cur_points] = x;
                DS[j].data[3*DS[j].cur_points + 1] = y;
                DS[j].data[3*DS[j].cur_points + 2] = z;
            }
            
            low = high;
            high += step;
        }
    }
}