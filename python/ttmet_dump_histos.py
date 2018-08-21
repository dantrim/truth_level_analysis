#!/usr/bin/env python

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
#r.gStyle.SetOptStat(False)
r.gROOT.SetBatch(True)

import os
import sys
import argparse


def get_inputs(args) :

    file_desc = args.input
    filenames = file_desc.split(",")
    return filenames

def dump_histograms(input_filename) :

    suffix = input_filename.split("_")[-1].replace(".root","")

    rfile = r.TFile.Open(input_filename)

    names = []
    for key in rfile.GetListOfKeys() :
        obj = key.ReadObj()
        if not obj.InheritsFrom("TH1") : continue
        names.append(obj.GetName())
        c = r.TCanvas("c_%s" % obj.GetName())
        c.cd()
        obj.SetLineColor(r.kBlack)
        obj.SetMinimum(0)
        obj.Draw("hist e")
        c.Update()
        c.SaveAs("./plots/histo_%s_%s.pdf" % (obj.GetName().replace("h_",""), suffix))

    return names

class RatioCanvas :
    def __init__(self, name) :
        self.canvas = r.TCanvas(name, name, 800, 600)
        self.upper_pad = r.TPad("upper", "upper", 0.0, 0.0, 1.0, 1.0)
        self.lower_pad = r.TPad("lower", "lower", 0.0, 0.0, 1.0, 1.0)
        self.set_pad_dims()

    def set_pad_dims(self) :

        can = self.canvas
        up  = self.upper_pad
        dn  = self.lower_pad

        can.cd()
        up_height = 0.75
        dn_height = 0.30
        up.SetPad(0.0, 1.0-up_height, 1.0, 1.0)
        dn.SetPad(0.0, 0.0, 1.0, dn_height)

        up.SetTickx(0)
        dn.SetGrid(0)
        dn.SetTicky(0)

        up.SetFrameFillColor(0)
        up.SetFillColor(0)

        # set right margins
        up.SetRightMargin(0.05)
        dn.SetRightMargin(0.05)

        # set left margins
        up.SetLeftMargin(0.14)
        dn.SetLeftMargin(0.14)

        # set top margins
        up.SetTopMargin(0.7 * up.GetTopMargin())
        
        # set bottom margins
        up.SetBottomMargin(0.09)
        dn.SetBottomMargin(0.4)

        up.Draw()
        dn.Draw()
        can.Update()

        self.canvas = can
        self.upper_pad = up
        self.lower_pad = dn

def dump_ratio_histograms(inputs, histo_names) :        

    r.gStyle.SetOptStat(False)

    suffixes = [ n.split("_")[-1].replace(".root","") for n in inputs ]
    colors = [r.kBlack, r.kRed]
    print "suffixes = %s" % suffixes

    histograms = {}
    for name in histo_names :
        histograms[name] = []

    rfiles = [r.TFile.Open(name) for name in inputs]
    for i, rfile in enumerate(rfiles) :
        for hname in histo_names :
            histo = rfile.Get(hname)
            histograms[hname].append(histo)

    for name in histo_names :
        title = name.replace("h_","")
        c = RatioCanvas(title)
        c.canvas.cd()
        c.upper_pad.cd()
        c.upper_pad.SetTicks(1,1)
        c.upper_pad.SetLogy(True)
        dont_draw = False

        val_to_scale = -1
        maxy = -1


        for ih, histo in enumerate(histograms[name]) :
            if not histo : break
            if histo.Integral() == 0 :
                dont_draw = True
                break
            histo.SetLineWidth(2)
            histo.SetLineColor(colors[ih])
            histo.Scale(1/histo.Integral())
            if ih == 0 :
                val_to_scale = histo.Integral()
            else :
                val_to_scale = histo.Integral() / val_to_scale
            histo.Scale(1/val_to_scale)


            #if histo.GetMaximum() > maxy :
            #    maxy = histo.GetMaximum()
            #histo.SetMaximum(0.4)

            yax = histo.GetYaxis()
            xax = histo.GetXaxis()

            yax.SetTitle("arb. units")
            yax.SetTitleOffset(1.5)
            xax.SetTitleOffset(100)
            xax.SetLabelOffset(100)
            option = "hist e"
            if ih != 0 : option += " same"

            #leg.AddEntry(histo, suffixes[ih], "l")
            histo.Draw(option)
            c.canvas.Update()

        if dont_draw :
            continue

        c.lower_pad.cd()
        c.lower_pad.SetTicks(1,1)

        if len(histograms[name]) != 2 : continue

        hratio = histograms[name][1].Clone("ratio_%s" % histograms[name][0].GetName())
        hratio.Divide(histograms[name][0])
        nbins = hratio.GetNbinsX()
        hratio.SetLineColor(colors[1])
        hratio.SetMinimum(0)
        hratio.SetMaximum(3)

        yax = hratio.GetYaxis()
        xax = hratio.GetXaxis()
        yax.SetTitle("Extras enabled / default")
        yax.SetTitleSize(2 * yax.GetTitleSize())
        yax.SetTitleOffset(0.7)
        yax.SetLabelSize(2 * yax.GetLabelSize())
        xax.SetTitleOffset(1.5)
        xax.SetTitleSize(3 * xax.GetTitleSize())
        xax.SetLabelSize(2*xax.GetLabelSize())
        xax.SetLabelOffset(0.05)

        hratio.Draw("hist e")
        #c.lower_pad.Update()
        c.canvas.Update()
            
        c.canvas.SaveAs("./plots/ratio_%s.pdf" % name)
            

def main() :

    parser = argparse.ArgumentParser(description = "Dumps plots from ttmet_truth_validation (or any other)")
    parser.add_argument("-i", "--input", required = True, \
        help = "Provide input file(s) (provide at most 2 comma separated files to make ratio plots")
    args = parser.parse_args()

    inputs = get_inputs(args)

    histo_names = []
    for input_filename in inputs :
        histo_names = dump_histograms(input_filename)
    if len(inputs) == 2 :
        dump_ratio_histograms(inputs, histo_names)


if __name__ == "__main__" :
    main()
