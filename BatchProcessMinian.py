# this is for batch processing of miniscope data
# 'conda activate minian' first in command line
# 'spyder' and then run this script; or you can just 'python BatchProcessMinian.py'
# MJH 02/27/2022

# 2.Setting up
## load modules
import itertools as itt
import os
import sys

import holoviews as hv
import numpy as np
import xarray as xr
from dask.distributed import Client, LocalCluster
from IPython.display import display
import zarr
import scipy.io as sio
import smtplib
from email.mime.text import MIMEText
## set path and parameters
animal = ['86']
state = ['Estrus']
sess = [['1','2'],['2','3'],['1','2'],['1','2'],['1','2'],['1','2'],['1','2']]
CNMFu2 = False
print('---Finish Setup---')
if __name__ == '__main__':
    for ii in range(np.size(animal)):
        for jj in range(np.size(state)):
            session = sess[ii]
            for kk in range(np.size(session)):
                dpath = 'H:/cxm/Estrus ds/AHN'+animal[ii]+'/'+state[jj]+'/Sess'+session[kk]+'/My_V4_Miniscope'
                SavePath = 'H:/cxm/Estrus ds/AHN'+animal[ii]+'/'+state[jj]+'/Sess'+session[kk]+'Res'
                intpath = dpath+'/minian_intermediate'
                if os.path.exists(dpath):
                    # Set up Initial Basic Parameters#
                    print(dpath)
                    minian_path = "."
                    minian_ds_path = os.path.join(dpath, "minian")
                    subset = dict(frame=slice(0, None))
                    subset_mc = None
                    interactive = False
                    output_size = 100
                    n_workers = int(os.getenv("MINIAN_NWORKERS", 6))
                    param_save_minian = {
                        "dpath": minian_ds_path,
                        "meta_dict": dict(session=-1, animal=-2),
                        "overwrite": True,
                    }
                    
                    # Pre-processing Parameters#
                    param_load_videos = {
                        "pattern": "[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]+\.avi$",
                        "dtype": np.uint8,
                        "downsample": dict(frame=1, height=1, width=1),
                        "downsample_strategy": "subset",
                    }
                    param_denoise = {"method": "median", "ksize": 7}
                    param_background_removal = {"method": "tophat", "wnd": 15}
                    
                    # Motion Correction Parameters#
                    subset_mc = None
                    param_estimate_motion = {"dim": "frame"}
                    
                    # Initialization Parameters#
                    param_seeds_init = {
                        "wnd_size": 1000,
                        "method": "rolling",
                        "stp_size": 500,
                        "max_wnd": 15,
                        "diff_thres": 3,
                    }
                    param_pnr_refine = {"noise_freq": 0.07, "thres": 1}
                    param_ks_refine = {"sig": 0.05}
                    param_seeds_merge = {"thres_dist": 10, "thres_corr": 0.8, "noise_freq": 0.07}
                    param_initialize = {"thres_corr": 0.8, "wnd": 10, "noise_freq": 0.07}
                    param_init_merge = {"thres_corr": 0.8}
                    
                    # CNMF Parameters#
                    param_get_noise = {"noise_range": (0.06, 0.5)}
                    param_first_spatial = {
                        "dl_wnd": 10,
                        "sparse_penal": 0.01,
                        "size_thres": (25, None),
                    }
                    param_first_temporal = {
                        "noise_freq": 0.07,
                        "sparse_penal": 1,
                        "p": 1,
                        "add_lag": 20,
                        "jac_thres": 0.2,
                    }
                    param_first_merge = {"thres_corr": 0.8}
                    param_second_spatial = {
                        "dl_wnd": 10,
                        "sparse_penal": 0.01,
                        "size_thres": (25, None),
                    }
                    param_second_temporal = {
                        "noise_freq": 0.07,
                        "sparse_penal": 1,
                        "p": 1,
                        "add_lag": 20,
                        "jac_thres": 0.4,
                    }
                    
                    os.environ["OMP_NUM_THREADS"] = "2"
                    os.environ["MKL_NUM_THREADS"] = "2"
                    os.environ["OPENBLAS_NUM_THREADS"] = "2"
                    os.environ["MINIAN_INTERMEDIATE"] = intpath
          ## import minian         
                    sys.path.append(minian_path)
                    from minian.cnmf import (
                        compute_AtC,
                        compute_trace,
                        get_noise_fft,
                        smooth_sig,
                        unit_merge,
                        update_spatial,
                        update_temporal,
                        update_background,
                    )
                    from minian.initialization import (
                        gmm_refine,
                        initA,
                        initC,
                        intensity_refine,
                        ks_refine,
                        pnr_refine,
                        seeds_init,
                        seeds_merge,
                    )
                    from minian.motion_correction import apply_transform, estimate_motion
                    from minian.preprocessing import denoise, remove_background
                    from minian.utilities import (
                        TaskAnnotation,
                        get_optimal_chk,
                        load_videos,
                        open_minian,
                        save_minian,
                    )
                    from minian.visualization import (
                        CNMFViewer,
                        VArrayViewer,
                        generate_videos,
                        visualize_gmm_fit,
                        visualize_motion,
                        visualize_preprocess,
                        visualize_seeds,
                        visualize_spatial_update,
                        visualize_temporal_update,
                        write_video,
                    )
                    print('---Module Minian Loaded---')
            ## module initialization
                    dpath = os.path.abspath(dpath)
                    hv.notebook_extension("bokeh", width=100)
            ## start cluster        
                    cluster = LocalCluster(
                        n_workers=n_workers,
                        memory_limit="16GB",
                        resources={"MEM": 1},
                        threads_per_worker=2,
                        dashboard_address=":8787",
                    )
                    annt_plugin = TaskAnnotation()
                    cluster.scheduler.add_plugin(annt_plugin)
                    client = Client(cluster)
# 3.Preprocessing
    ## loading videos and visualization                    
                    param_load_videos
                    
                    varr = load_videos(dpath, **param_load_videos)
                    chk, _ = get_optimal_chk(varr, dtype=float)
                    
                    varr = save_minian(
                        varr.chunk({"frame": chk["frame"], "height": -1, "width": -1}).rename("varr"),
                        intpath,
                        overwrite=True,
                    )
    
                    # varr.sel(frame=slice(0,799)) #subset part of video
                    varr_ref = varr.sel(subset)
                    varr_min = varr_ref.min("frame").compute()
                    varr_ref = varr_ref - varr_min
                    
                    print('---Now Performing Denoise---')
                    print(dpath,'denoising, may take a while')
                    param_denoise
                    varr_ref = denoise(varr_ref, **param_denoise)
                    print('---Now Removing Background---')
                    param_background_removal
                    varr_ref = remove_background(varr_ref, **param_background_removal)
                    varr_ref = save_minian(varr_ref.rename("varr_ref"), dpath=intpath, overwrite=True)
                    
                    param_estimate_motion
                    motion = estimate_motion(varr_ref.sel(subset_mc), **param_estimate_motion)
                    param_save_minian
                    print('---Now Saving after Removing Background---')
                    motion = save_minian(
                        motion.rename("motion").chunk({"frame": chk["frame"]}), **param_save_minian
                    )
                    Y = apply_transform(varr_ref, motion, fill=0)
                    Y_fm_chk = save_minian(Y.astype(float).rename("Y_fm_chk"), intpath, overwrite=True)
                    Y_hw_chk = save_minian(
                        Y_fm_chk.rename("Y_hw_chk"),
                        intpath,
                        overwrite=True,
                        chunks={"frame": -1, "height": chk["height"], "width": chk["width"]},
                    )
                    print('---Now Setting Dict---')
                    # im_opts = dict(
                    #     frame_width=500,
                    #     aspect=varr_ref.sizes["width"] / varr_ref.sizes["height"],
                    #     cmap="Viridis",
                    #     colorbar=True,
                    # )
                    # (
                    #     regrid(
                    #         hv.Image(
                    #             varr_ref.max("frame").compute().astype(np.float32),
                    #             ["width", "height"],
                    #             label="before_mc",
                    #         ).opts(**im_opts)
                    #     )
                    #     + regrid(
                    #         hv.Image(
                    #             Y_hw_chk.max("frame").compute().astype(np.float32),
                    #             ["width", "height"],
                    #             label="after_mc",
                    #         ).opts(**im_opts)
                    #     )
                    # )
                    # vid_arr = xr.concat([varr_ref, Y_fm_chk], "width").chunk({"width": -1})
                    # print('---Writing Motion Correction Video Now---')
                    # write_video(vid_arr, "minian_mc.mp4", dpath)
                    
                    max_proj = save_minian(
                        Y_fm_chk.max("frame").rename("max_proj"), **param_save_minian
                    ).compute()
                    
                    print('---SeedsInit---')
                    param_seeds_init
                    seeds = seeds_init(Y_fm_chk, **param_seeds_init)
                    seeds.head() 
                    if interactive:
                        noise_freq_list = [0.005, 0.01, 0.02, 0.06, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8]
                        example_seeds = seeds.sample(6, axis="rows")
                        example_trace = Y_hw_chk.sel(
                            height=example_seeds["height"].to_xarray(),
                            width=example_seeds["width"].to_xarray(),
                        ).rename(**{"index": "seed"})
                        smooth_dict = dict()
                        for freq in noise_freq_list:
                            trace_smth_low = smooth_sig(example_trace, freq)
                            trace_smth_high = smooth_sig(example_trace, freq, btype="high")
                            trace_smth_low = trace_smth_low.compute()
                            trace_smth_high = trace_smth_high.compute()
                            hv_trace = hv.HoloMap(
                                {
                                    "signal": (
                                        hv.Dataset(trace_smth_low)
                                        .to(hv.Curve, kdims=["frame"])
                                        .opts(frame_width=300, aspect=2, ylabel="Signal (A.U.)")
                                    ),
                                    "noise": (
                                        hv.Dataset(trace_smth_high)
                                        .to(hv.Curve, kdims=["frame"])
                                        .opts(frame_width=300, aspect=2, ylabel="Signal (A.U.)")
                                    ),
                                },
                                kdims="trace",
                            ).collate()
                            smooth_dict[freq] = hv_trace
                            
                    param_pnr_refine
                    seeds, pnr, gmm = pnr_refine(Y_hw_chk, seeds, **param_pnr_refine)
                    seeds.head()
                    if gmm:
                        display(visualize_gmm_fit(pnr, gmm, 100))
                    else:
                        print("nothing to show")
                    
                    param_ks_refine
                    seeds = ks_refine(Y_hw_chk, seeds, **param_ks_refine)
                    param_seeds_merge
                    seeds_final = seeds[seeds["mask_ks"] & seeds["mask_pnr"]].reset_index(drop=True)
                    seeds_final = seeds_merge(Y_hw_chk, max_proj, seeds_final, **param_seeds_merge)
                    
                    print('---Now Doing Param_initialize---')
                    param_initialize
                    A_init = initA(Y_hw_chk, seeds_final[seeds_final["mask_mrg"]], **param_initialize)
                    A_init = save_minian(A_init.rename("A_init"), intpath, overwrite=True)
                    C_init = initC(Y_fm_chk, A_init)
                    C_init = save_minian(
                        C_init.rename("C_init"), intpath, overwrite=True, chunks={"unit_id": 1, "frame": -1}
                    )
                    
                    param_init_merge
                    try:
                        A, C = unit_merge(A_init, C_init, **param_init_merge)
                    except ValueError as e1:
                        A = A_init
                        C = C_init
                    A = save_minian(A.rename("A"), intpath, overwrite=True)
                    C = save_minian(C.rename("C"), intpath, overwrite=True)
                    C_chk = save_minian(
                        C.rename("C_chk"),
                        intpath,
                        overwrite=True,
                        chunks={"unit_id": -1, "frame": chk["frame"]},
                    )
                    b, f = update_background(Y_fm_chk, A, C_chk)
                    f = save_minian(f.rename("f"), intpath, overwrite=True)
                    b = save_minian(b.rename("b"), intpath, overwrite=True)
                    
                    # CNMF
                    print('--- CNMF running ---')
                    print('--- Now Working on',dpath,'---')
                    param_get_noise
                    sn_spatial = get_noise_fft(Y_hw_chk, **param_get_noise)
                    sn_spatial = save_minian(sn_spatial.rename("sn_spatial"), intpath, overwrite=True)
                    if interactive:
                        units = np.random.choice(A.coords["unit_id"], 10, replace=False)
                        units.sort()
                        A_sub = A.sel(unit_id=units).persist()
                        C_sub = C.sel(unit_id=units).persist()
                    if interactive:
                        sprs_ls = [0.005, 0.01, 0.05]
                        A_dict = dict()
                        C_dict = dict()
                        for cur_sprs in sprs_ls:
                            cur_A, cur_mask, cur_norm = update_spatial(
                                Y_hw_chk,
                                A_sub,
                                C_sub,
                                sn_spatial,
                                in_memory=True,
                                dl_wnd=param_first_spatial["dl_wnd"],
                                sparse_penal=cur_sprs,
                            )
                            if cur_A.sizes["unit_id"]:
                                A_dict[cur_sprs] = cur_A.compute()
                                C_dict[cur_sprs] = C_sub.sel(unit_id=cur_mask).compute()
                        hv_res = visualize_spatial_update(A_dict, C_dict, kdims=["sparse penalty"])
                    
                    param_first_spatial
                    A_new, mask, norm_fac = update_spatial(
                        Y_hw_chk, A, C, sn_spatial, **param_first_spatial
                    )
                    C_new = save_minian(
                        (C.sel(unit_id=mask) * norm_fac).rename("C_new"), intpath, overwrite=True
                    )
                    C_chk_new = save_minian(
                        (C_chk.sel(unit_id=mask) * norm_fac).rename("C_chk_new"), intpath, overwrite=True
                    )
                    b_new, f_new = update_background(Y_fm_chk, A_new, C_chk_new)
                    #save1
                    print('---Saving After first Spatial Update---')
                    A = save_minian(
                        A_new.rename("A"),
                        intpath,
                        overwrite=True,
                        chunks={"unit_id": 1, "height": -1, "width": -1},
                    )
                    b = save_minian(b_new.rename("b"), intpath, overwrite=True)
                    f = save_minian(
                        f_new.chunk({"frame": chk["frame"]}).rename("f"), intpath, overwrite=True
                    )
                    C = save_minian(C_new.rename("C"), intpath, overwrite=True)
                    C_chk = save_minian(C_chk_new.rename("C_chk"), intpath, overwrite=True)
                    #
                    if interactive:
                        units = np.random.choice(A.coords["unit_id"], 10, replace=False)
                        units.sort()
                        A_sub = A.sel(unit_id=units).persist()
                        C_sub = C_chk.sel(unit_id=units).persist()
                    if interactive:
                        p_ls = [1]
                        sprs_ls = [0.1, 0.5, 1, 2]
                        add_ls = [20]
                        noise_ls = [0.06]
                        YA_dict, C_dict, S_dict, g_dict, sig_dict, A_dict = [dict() for _ in range(6)]
                        YrA = (
                            compute_trace(Y_fm_chk, A_sub, b, C_sub, f)
                            .persist()
                            .chunk({"unit_id": 1, "frame": -1})
                        )
                        for cur_p, cur_sprs, cur_add, cur_noise in itt.product(
                            p_ls, sprs_ls, add_ls, noise_ls
                        ):
                            ks = (cur_p, cur_sprs, cur_add, cur_noise)
                            print(
                                "p:{}, sparse penalty:{}, additional lag:{}, noise frequency:{}".format(
                                    cur_p, cur_sprs, cur_add, cur_noise
                                )
                            )
                            cur_C, cur_S, cur_b0, cur_c0, cur_g, cur_mask = update_temporal(
                                A_sub,
                                C_sub,
                                YrA=YrA,
                                sparse_penal=cur_sprs,
                                p=cur_p,
                                use_smooth=True,
                                add_lag=cur_add,
                                noise_freq=cur_noise,
                            )
                            YA_dict[ks], C_dict[ks], S_dict[ks], g_dict[ks], sig_dict[ks], A_dict[ks] = (
                                YrA.compute(),
                                cur_C.compute(),
                                cur_S.compute(),
                                cur_g.compute(),
                                (cur_C + cur_b0 + cur_c0).compute(),
                                A_sub.compute(),
                            )
                        hv_res = visualize_temporal_update(
                            YA_dict,
                            C_dict,
                            S_dict,
                            g_dict,
                            sig_dict,
                            A_dict,
                            kdims=["p", "sparse penalty", "additional lag", "noise frequency"],
                        )
                    YrA = save_minian(
                        compute_trace(Y_fm_chk, A, b, C_chk, f).rename("YrA"),
                        intpath,
                        overwrite=True,
                        chunks={"unit_id": 1, "frame": -1},
                    )
                    
                    param_first_temporal
                    C_new, S_new, b0_new, c0_new, g, mask = update_temporal(
                        A, C, YrA=YrA, **param_first_temporal
                    )
                    #save2
                    print('---Saving After First Temproal Update---')
                    C = save_minian(
                        C_new.rename("C").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
                    )
                    C_chk = save_minian(
                        C.rename("C_chk"),
                        intpath,
                        overwrite=True,
                        chunks={"unit_id": -1, "frame": chk["frame"]},
                    )
                    S = save_minian(
                        S_new.rename("S").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
                    )
                    b0 = save_minian(
                        b0_new.rename("b0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
                    )
                    c0 = save_minian(
                        c0_new.rename("c0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
                    )
                    A = A.sel(unit_id=C.coords["unit_id"].values)
                    
                    print('---Merge Units---')
                    print('Now Processing: ',dpath,)
                    param_first_merge
                    try:
                        A_mrg, C_mrg, [sig_mrg] = unit_merge(A, C, [C + b0 + c0], **param_first_merge)
                    except ValueError as e:
                        A_mrg = A
                        C_mrg = C
                        sig_mrg = sig = C + b0 + c0
                        print('No merge Now')
                    A = save_minian(A_mrg.rename("A_mrg"), intpath, overwrite=True)
                    C = save_minian(C_mrg.rename("C_mrg"), intpath, overwrite=True)
                    C_chk = save_minian(
                        C.rename("C_mrg_chk"),
                        intpath,
                        overwrite=True,
                        chunks={"unit_id": -1, "frame": chk["frame"]},
                    )
                    sig = save_minian(sig_mrg.rename("sig_mrg"), intpath, overwrite=True)
                    
                    print('---First CNMF Update Ends Here---')
                    if CNMFu2:
   
                        A_new, mask, norm_fac = update_spatial(
                            Y_hw_chk, A, C, sn_spatial, **param_second_spatial
                        )
                        C_new = save_minian(
                            (C.sel(unit_id=mask) * norm_fac).rename("C_new"), intpath, overwrite=True
                        )
                        C_chk_new = save_minian(
                            (C_chk.sel(unit_id=mask) * norm_fac).rename("C_chk_new"), intpath, overwrite=True
                        )
                        
                        b_new, f_new = update_background(Y_fm_chk, A_new, C_chk_new)
                        A = save_minian(
                            A_new.rename("A"),
                            intpath,
                            overwrite=True,
                            chunks={"unit_id": 1, "height": -1, "width": -1},
                        )
                        b = save_minian(b_new.rename("b"), intpath, overwrite=True)
                        f = save_minian(
                            f_new.chunk({"frame": chk["frame"]}).rename("f"), intpath, overwrite=True
                        )
                        C = save_minian(C_new.rename("C"), intpath, overwrite=True)
                        C_chk = save_minian(C_chk_new.rename("C_chk"), intpath, overwrite=True)
                        
                        
                        YrA = save_minian(
                            compute_trace(Y_fm_chk, A, b, C_chk, f).rename("YrA"),
                            intpath,
                            overwrite=True,
                            chunks={"unit_id": 1, "frame": -1},
                        )
                        C_new, S_new, b0_new, c0_new, g, mask = update_temporal(
                            A, C, YrA=YrA, **param_second_temporal
                        )


# 6.3.6 save result
                        C = save_minian(
                            C_new.rename("C").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
                        )
                        C_chk = save_minian(
                            C.rename("C_chk"),
                            intpath,
                            overwrite=True,
                            chunks={"unit_id": -1, "frame": chk["frame"]},
                        )
                        S = save_minian(
                            S_new.rename("S").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
                        )
                        b0 = save_minian(
                            b0_new.rename("b0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
                        )
                        c0 = save_minian(
                            c0_new.rename("c0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
                        )
                        A = A.sel(unit_id=C.coords["unit_id"].values)
                    print('---Now Generating Final Video, Saving---')
                    generate_videos(varr.sel(subset), Y_fm_chk, A=A, C=C_chk, vpath=dpath)

                    #final save
                    print('---Second CNMF Ends Here, Saving---')
                    A = save_minian(A.rename("A"), **param_save_minian)
                    C = save_minian(C.rename("C"), **param_save_minian)
                    S = save_minian(S.rename("S"), **param_save_minian)
                    c0 = save_minian(c0.rename("c0"), **param_save_minian)
                    b0 = save_minian(b0.rename("b0"), **param_save_minian)
                    b = save_minian(b.rename("b"), **param_save_minian)
                    f = save_minian(f.rename("f"), **param_save_minian)
# 7. zarr to matlab                    
                    def readzarr(filename,outmat,ani,SavePath):
                        zarray = zarr.open(filename + outmat + '.zarr')
                        out = zarray[outmat]
                        outarray = out[:]
                        if os.path.exists(SavePath):
                            print('Path Made')
                        else:
                            os.mkdir(SavePath)
                        sio.savemat(SavePath+'/'+ani+outmat+'.mat',{'array':outarray})
                        return
                    dspMat = SavePath
                    AniName = animal[ii]
                    ZarrFile = dpath+'/minian/'
                    readzarr(ZarrFile,'C',AniName,dspMat)
                    readzarr(ZarrFile,'A',AniName,dspMat)
                    readzarr(ZarrFile,'b0',AniName,dspMat)
                    readzarr(ZarrFile,'c0',AniName,dspMat)
                    readzarr(ZarrFile,'max_proj',AniName,dspMat)

                    Varray = zarr.open(intpath+'/Y_fm_chk.zarr')
                    out = Varray['Y_fm_chk']
                    outarray = out.astype(np.uint8)
                    Shape = np.shape(outarray[:])
                    length = Shape[0]
                    split_mat = np.arange(0,length,9000)
                    if length>9000:
                        for idx in range(np.size(split_mat)-1):
                            print('Now Working on',idx,'th mat')
                            start = split_mat[idx]
                            end = split_mat[idx+1] - 1
                            sio.savemat(SavePath+'/varr'+str(idx+1)+'.mat',{'array':outarray[start:end,:,:]})
                        sio.savemat(SavePath+'/varr'+str(idx+2)+'.mat',{'array':outarray[split_mat[idx+1]:length,:,:]})    
                        
                        print('----',dpath,'Done! ----')
                    else:
                        sio.savemat(SavePath+'/varr'+str(1)+'.mat',{'array':outarray[:,:,:]})
                        
                    client.close()
                    cluster.close()
                    mail_host = 'smtp.qq.com'  
                    mail_user = '1132073756'
                    mail_pass = 'zzyfapaqqrdnjhde'   

                    sender = '1132073756@qq.com'  
                    receivers = ['1132073756@qq.com']  

                    message = MIMEText('Come Here To Check','plain','utf-8')     
                    message['Subject'] = 'Your Program Has Been Done' 
                    message['From'] = sender     
                    message['To'] = receivers[0]  

                    try:
                        smtpObj = smtplib.SMTP() 
                        smtpObj.connect(mail_host,25)
                        smtpObj.login(mail_user,mail_pass) 
                        smtpObj.sendmail(
                            sender,receivers,message.as_string()) 
                        smtpObj.quit() 
                        print('Done! Successfully Sending the Email')
                    except smtplib.SMTPException as e:
                        print('error',e)
                else:
                    print('---CANNOT FIND YOUR MINISCOPE DATA AT YOUR LOCATION, PLEASE CHECK!!!---')
                    continue




