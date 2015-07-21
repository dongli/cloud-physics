load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    root = "."
    start_time = 0
    end_time = 1000
    time_step = 10
    time_stamp = "%5.5d"

    file_prefix = "output"

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '"+root+"/"+file_prefix+"."+time_stamp+".nc\n' $i; done")

    wks = gsn_open_wks("pdf", "droplet_spectrum")

    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@xyLineThicknessF = 3
    res@trYMaxF = 100
    res@trXMinF = 0.00000001
    res@trXMaxF = 100
    res@trXLog = True
    res@tiXAxisString = "Droplet Radius (cm)"
    res@tiYAxisString = "Drop-size Distribution (cm-4)"
    res@tiMainString = "      Condensation Growth ~C~Solved by Lagrangian Method"

    res_text = True

    do i = 0, dimsizes(fs)-1
        f = addfile(fs(i), "r")
        if (.not. isMonotonic(f->r)) then
            system("echo '[Error]: r is not monotonic!'")
        end if
        system("echo 'plot file "+fs(i)+" ...'")
        plot = gsn_csm_xy(wks, f->r, f->f, res)
        getvalues plot
            "tmYLLabelFontHeightF"   : fheight
            "tmXTValues"             : tmXTValues
            "tmYLValues"             : tmYLValues
        end getvalues
        nx = dimsizes(tmXTValues)
        ny = dimsizes(tmYLValues)
        res_text@txFontHeightF = fheight
        text = gsn_add_text(wks, plot, "~F8~D~F0~t = 0.01s", 2.0e-8*tmXTValues(nx-1), 0.85*tmYLValues(ny-1), res_text)
        draw(plot)
        frame(wks)
    end do

end
