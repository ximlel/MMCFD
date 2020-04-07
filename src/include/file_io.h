/* This function reads the initial data file. The function 
 * initialize return a pointer pointing to the position of
 * a block of memory consisting (m+1) variables* of type
 * double. The value of first of these variables is m.
 * The following m variables are the initial value.
 */
int flu_var_init(char * add, double * F, const int r_or_c);
/* This function read the configuration data file,
 * and store the configuration data in the array "config".
 * config[?] stand for what? (see  Configuration_instructions.pdf)
 */
int configurate(char * add);
void example_io(const char *example, char *add_mkdir, const int i_or_o);
struct flu_var flu_conf_load(const char *example);


void file_write_TEC(const struct flu_var FV, const struct mesh_var mv, const char * problem, const double time, const int dim);
void file_write_VTK_3D(const struct flu_var FV, const struct mesh_var mv, const char * problem);
