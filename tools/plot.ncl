load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    root = "."
    start_time = 0
    end_time = 30000
    time_step = 300
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
    res@xyMarkLineMode = "MarkLines"
    res@trYMaxF = 100
    res@trXMinF = 0
    res@trXMaxF = 50
    res@tiXAxisString = "Droplet Radius (~F8~m~F25~m)"
    res@tiYAxisString = "Drop-size Distribution (cm~S~-4~N~)"
    res@tiMainString = "      Condensation Growth ~C~Solved by Lagrangian Method"

    res_text = True

    do i = 0, dimsizes(fs)-1
        f = addfile(fs(i), "r")
        if (.not. isMonotonic(f->r)) then
            system("echo '[Error]: r is not monotonic!'")
            exit
        end if
        system("echo 'plot file "+fs(i)+" ...'")
        plot = gsn_csm_xy(wks, f->r*1.0d4, f->f, res)
        getvalues plot
            "tmYLLabelFontHeightF"   : fheight
            "tmXTValues"             : tmXTValues
            "tmYLValues"             : tmYLValues
        end getvalues
        nx = dimsizes(tmXTValues)
        ny = dimsizes(tmYLValues)
        res_text@txFontHeightF = fheight
        text = gsn_add_text(wks, plot, "~F8~D~F25~t = "+sprintf("%4.1f", f@time_step_size)+"s", 0.2*tmXTValues(nx-1), 0.85*tmYLValues(ny-1), res_text)
        text = gsn_add_text(wks, plot, "#bin = "+dimsizes(f->r), 0.2*tmXTValues(nx-1), 0.78*tmYLValues(ny-1), res_text)
        text = gsn_add_text(wks, plot, "f~B~max~N~ = "+sprintf("%3.0f", max(f->f)), 0.2*tmXTValues(nx-1), 0.71*tmYLValues(ny-1), res_text)
        draw(plot)
        frame(wks)
        delete(tmXTValues)
        delete(tmYLValues)
        print(min(f->dr))
    end do

end
