// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos
// VERSION 2022/09/01

#include <unistd.h>
#include <groan.h>

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;

/*
 * Parses command line arguments.
 * Returns zero, if parsing has been successful. Else returns non-zero.
 */
int get_arguments(
        int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **output_file,
        char **reference_atoms,
        int *skip,
        int *center_x,
        int *center_y,
        int *center_z) 
{
    int gro_specified = 0, output_specified = 0;

    int opt = 0;
    while((opt = getopt(argc, argv, "c:f:n:o:r:s:xyzh")) != -1) {
        switch (opt) {
        // help
        case 'h':
            return 1;
        // gro file to read
        case 'c':
            *gro_file = optarg;
            gro_specified = 1;
            break;
        // xtc file to read
        case 'f':
            *xtc_file = optarg;
            break;
        // ndx file to read
        case 'n':
            *ndx_file = optarg;
            break;
        // output file name
        case 'o':
            *output_file = optarg;
            output_specified = 1;
            break;
        // reference atoms
        case 'r':
            *reference_atoms = optarg;
            break;
        // skip frames
        case 's':
            if (sscanf(optarg, "%d", skip) != 1) {
                fprintf(stderr, "Could not parse skip value (flag '-s').\n");
                return 1;
            }
            
            if (*skip <= 0) {
                fprintf(stderr, "Skip must be positive.\n");
                return 1;
            }
            break;
        // centering in individual dimensions
        case 'x':
            *center_x = 1;
            break;
        case 'y':
            *center_y = 1;
            break;
        case 'z':
            *center_z = 1;
            break;
        default:
            //fprintf(stderr, "Unknown command line option: %c.\n", opt);
            return 1;
        }
    }

    if (!gro_specified || !output_specified) {
        fprintf(stderr, "Gro file and output file must always be supplied.\n");
        return 1;
    }
    return 0;
}

void print_usage(const char *program_name)
{
    printf("Usage: %s -c GRO_FILE -o OUTPUT_FILE [OPTION]...\n", program_name);
    printf("\nOPTIONS\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read (optional)\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-o STRING        output file name\n");
    printf("-r STRING        selection of atoms centered (default: Protein)\n");
    printf("-s INTEGER       only center every Nth frame (default: 1)\n");
    printf("-x/-y/-z         center in individual x/y/z dimensions (default: center in xyz)\n");
    printf("\n");
}

inline void set_translation(vec_t translation, box_t box, vec_t center, const int x, const int y, const int z)
{
    if (x) translation[0] = (box[0] / 2) - center[0];
    if (y) translation[1] = (box[1] / 2) - center[1];
    if (z) translation[2] = (box[2] / 2) - center[2];
}

int main(int argc, char **argv)
{
    // get arguments
    char *gro_file = NULL;
    char *xtc_file = NULL;
    char *ndx_file = "index.ndx";
    char *output_file = NULL;
    char *reference_atoms = "Protein";
    int skip = 1;
    int center_x = 0;
    int center_y = 0;
    int center_z = 0;

    if (get_arguments(argc, argv, &gro_file, &xtc_file, &ndx_file, &output_file, &reference_atoms, &skip, &center_x, &center_y, &center_z) != 0) {
        print_usage(argv[0]);
        return 1;
    }

    // check that the paths to input and output files are not the same
    // this does not work if the paths are different but point to the same file!
    if (!strcmp(gro_file, output_file)) {
        fprintf(stderr, "Input gro file %s and output file %s are the same file.\n", gro_file, output_file);
        return 1;
    }

    if (xtc_file != NULL) {
        if (!strcmp(gro_file, xtc_file)) {
            fprintf(stderr, "Input gro file %s and input xtc file %s are the same file.\n", gro_file, xtc_file);
            return 1;
        }

        if (!strcmp(xtc_file, output_file)) {
            fprintf(stderr, "Input xtc file %s and output file %s are the same file.\n", xtc_file, output_file);
            return 1;
        }
    }

    // if no center dimension has been selected, use all of them
    if (!center_x && !center_y && !center_z) {
        center_x = 1;
        center_y = 1;
        center_z = 1;
    }

    // read gro file
    system_t *system = load_gro(gro_file);
    if (system == NULL) return 1;

    // try reading ndx file (ignore if this fails)
    dict_t *ndx_groups = read_ndx(ndx_file, system);

    // select all atoms
    atom_selection_t *all = select_system(system);

    // select reference atoms
    select_t *reference = smart_select(all, reference_atoms, ndx_groups);
    if (reference == NULL || reference->n_atoms == 0) {
        fprintf(stderr, "No reference atoms ('%s') found.\n", reference_atoms);

        dict_destroy(ndx_groups);
        free(system);
        free(all);
        free(reference);
        return 1;
    }

    // if there is no xtc file supplied, just center gro file and write it
    if (xtc_file == NULL) {
        FILE *output = fopen(output_file, "w");
        if (output == NULL) {
            fprintf(stderr, "File %s could not be opened for writing.\n", output_file);
            dict_destroy(ndx_groups);
            free(system);
            free(all);
            free(reference);
            return 1;
        }

        vec_t center = {0.0f};

        center_of_geometry(reference, center, system->box);
        vec_t translation = {0.0f};
        set_translation(translation, system->box, center, center_x, center_y, center_z);
        selection_translate(all, translation, system->box);

        int return_code = 0;
        if (write_gro(output, all, system->box, velocities, "Generated using `center`.") != 0) {
            fprintf(stderr, "Writing has failed.\n");
            return_code = 1;
        }
        
        dict_destroy(ndx_groups);
        free(system);
        free(all);
        free(reference);
        fclose(output);
        return return_code;
    }

    // open xtc file for reading
    XDRFILE *xtc = xdrfile_open(xtc_file, "r");
    if (xtc == NULL) {
        fprintf(stderr, "File %s could not be read as an xtc file.\n", xtc_file);
        dict_destroy(ndx_groups);
        free(system);
        free(all);
        free(reference);
        return 1;
    }

    // set all velocities of all particles to zero (xtc does not contain velocities)
    reset_velocities(system);

    // check that the gro file and the xtc file match each other
    if (!validate_xtc(xtc_file, (int) system->n_atoms)) {
        fprintf(stderr, "Number of atoms in %s does not match %s.\n", xtc_file, gro_file);
        xdrfile_close(xtc);
        dict_destroy(ndx_groups);
        free(system);
        free(all);
        free(reference);
        return 1;
    }

    // open output xtc for writing
    XDRFILE *output = xdrfile_open(output_file, "w");
    if (output == NULL) {
        fprintf(stderr, "File %s could not be opened for writing.\n", output_file);
        xdrfile_close(xtc);
        dict_destroy(ndx_groups);
        free(system);
        free(all);
        free(reference);
        return 1;
    }

    // loop through input xtc file, center each frame and write it into output
    int frames = 0;
    while (read_xtc_step(xtc, system) == 0) {

        // print info about the progress of reading and writing
        if ((int) system->time % PROGRESS_FREQ == 0) {
            printf("Step: %d. Time: %.0f ps\r", system->step, system->time);
            fflush(stdout);
        }

        if (frames % skip != 0) {
            ++frames;
            continue;
        }

        vec_t center = {0.0f};

        center_of_geometry(reference, center, system->box);
        vec_t translation = {0.0f};
        set_translation(translation, system->box, center, center_x, center_y, center_z);
        selection_translate(all, translation, system->box);

        if (write_xtc_step(output, all, system->step, system->time, system->box, system->precision) == 1) {
            fprintf(stderr, "Writing has failed.\n");
            dict_destroy(ndx_groups);
            free(all);
            free(reference);
            xdrfile_close(xtc);
            xdrfile_close(output);
            return 1;
        }

        ++frames;
    }
    printf("\n");

    dict_destroy(ndx_groups);
    free(all);
    free(reference);

    xdrfile_close(xtc);
    xdrfile_close(output);

    free(system);
    return 0;
}