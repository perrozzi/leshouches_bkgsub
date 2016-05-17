#!/bin/sh 
rivet-mkhtml -c style.plot -m/WWbb/WBF_njets.. WWbbLO.yoda:Scale=677.2:Title="WWbb" ttLO.yoda:Scale=416.3:Title="tt" tWb.yoda:Scale=7.481 WWbb_notopLO.yoda:Scale=1.168:Title="WWbb (no top)" 
cp -r plots ~/public/html/bkgsub_sub/
rivet-mkhtml -c style.plot -m/WWbb/WBF_njets.. WWbbLO.yoda:Scale=677.2:Title="WWbb" tt-WWbb_notop-tWb.yoda:Title="WWbb (no interference)"
cp -r plots ~/public/html/bkgsub_sum/
