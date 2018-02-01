void thresholding(CollocationPoint** collPnt, double epsilon);
void adjacent_zone(CollocationPoint** collPnt,int buffer_width,int buffer_height);
void reconstruction_check(CollocationPoint** collPnt);
void compress(CollocationPoint** collPnt, double** gridpts, double** funcpts, int** spatial_identifiers, int** level_identifiers, int& counter);
void decompress(CollocationPoint** collPnt, double* funcpts, int* spatial_identifier, int* level_identifier, int nactive);
