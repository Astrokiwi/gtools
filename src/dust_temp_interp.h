struct coolHeatDust {
    double dCool;
    double dHeat;
    double dustT;
    double opac_abs;
    double opac_scat;
    double dg;
    double column_out;
    double arad;
    
    int id,it,ii,ic;
    double fd,ft,fi,fc;
};


class CoolHeatTab {
    public:
        struct coolHeatDust interpTab(double density, double temperature, double intensity, double column_in);

        //CoolHeatTab(std::string flabels,std::string ftab);
        CoolHeatTab(const char* flabels,const char* ftab);

    private:
        int agn_tab_index(int id, int it, int ii, int is);
        int value_to_index(double edges[], int nedges, double value);
        double right_weight(double edges[], int nedges, int index, double value);
        
        // heating/cooling tables under AGN (table values)
        double  *agn_heat_tab,
                *agn_cool_tab,
                *agn_dust_tab,
                *agn_dg_tab,
                *agn_opac_scat_tab,
                *agn_opac_abs_tab,
                *agn_arad_tab,
                *agn_column_out_tab;
        double *tables[4];
        int *ntabs[4];


        // table axes
        int agn_ntemp,agn_ndense,agn_nintensity,agn_ncolumn_in;
        double *agn_temp_vals,*agn_dense_vals,*agn_intensity_vals,*agn_column_in_vals;
        
        bool data_loaded;
        
};
