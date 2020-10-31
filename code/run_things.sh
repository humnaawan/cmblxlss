#!/bin/bash

grecon=0 # if running with fgs

if [ ${grecon} == 0 ] ;
then
    for cmb_data_type in 'mv' 'TT'
    do
        python save_correlated_alms.py \
            --cmb_data_type=${cmb_data_type}

        python get_spectra.py \
            --cmb_data_type=${cmb_data_type}
    done
fi

if [ ${grecon} == 1 ] ;
then
    for cmb_data_type in 'mv' 'TT'
    do
        python save_correlated_alms.py \
            --grecon --cmb_data_type=${cmb_data_type}

        python get_spectra.py \
            --grecon --cmb_data_type=${cmb_data_type}
    done
fi





