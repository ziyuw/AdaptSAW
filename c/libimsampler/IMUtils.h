void initIMState(IMState* p_im_state, const IMInputModel* p_in_model,
                 const IMSimParams* p_sim_params);
void freeIMState(IMState* p_im_state, const IMInputModel* p_in_model,
                 const IMSimParams* p_sim_params);
void setIMState(IMState* p_im_state,
                const IMInputModel* p_in_model,
                const int *S_new, const double* h_effective_new);

void initSAWData(SAWData* p_saw_data, const IMInputModel* p_in_model,
                 const IMSimParams* p_sim_params);
void freeSAWData(SAWData* p_saw_data, const IMInputModel* p_in_model,
                 const IMSimParams* p_sim_params);

void proposeIMMove(SAWData* p_saw_data,
                   IMState* p_im_state, const IMInputModel* p_in_model,
                   const IMSimParams* p_sim_params);
void undoIMMove(SAWData* p_saw_data,
                IMState* p_im_state, const IMInputModel* p_in_model,
                const IMSimParams* p_sim_params);

double compute_E_in(const int* S_in,
                    const double* h_effective_in, const double* h,
                    int num_nodes);

void free_input_model(IMInputModel *in_model);

