void thresholding(CollocationPoint** collPnt, double epsilon);
void reconstruction(CollocationPoint** collPnt);
void compute_derivatives( CollocationPoint** collPnt ); 
void time_integrate(CollocationPoint** collPnt,double h,double c,double alpha);
void seed_grid(CollocationPoint** collPnt);
void extend_mask(CollocationPoint** collPnt,int buffer_width);
