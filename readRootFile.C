#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>



double Tetsuo_tuning(double *Radius, double *beta){
        const double hbarc = 197.3269602;
        const double  coeff[9] = {0.14, -13.76, 306.81, -2729.49, 11744.48, -26288.42, 30558.08, -16730.80, 3051.85};
        const double coeffc[9] = {0.07, -6.48, 131.39, -1021.60, 3548.93, -5159.25, 667.40, 5175.64, -3446.89};
        double coeff_rest[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
        // here we sustract the entries of the two vectors coeff and coeffc
        double mn[9] = {100., 200., 300.,400.,500.,600.,700.,800.,900.};
        for(int i = 0 ; i < 9 ; i++){
        mn[i] = mn[i]/hbarc;    
}

        for( int j = 0 ; j < 9 ; j++){
        coeff_rest[j] = (coeff[j] - coeffc[j]);
}
        const double pi = 3.1415926535897932384626433832795028841971693993751;
        double Result = 0;
        const double lambda = 1000./hbarc;
        const double r=Radius[0];
        double beta1 = beta[0];
        double V_rest = 0;
        double V_c = 0;
        //const double& gamma = Radius[3];
//      printf("\n gamma is: %.2f", gamma);
        for (unsigned uPar = 0; uPar < 9; uPar++){
        V_rest += hbarc*(((1./(4.*pi*r))*coeff_rest[uPar]*(pow(pow(lambda,2) /
         (pow(lambda,2) - pow(mn[uPar],2)),2)))*(exp(-mn[uPar]*r)-exp(-lambda*r)*
        ((pow(lambda,2) - pow(mn[uPar],2))*r + 2*lambda) / (2*lambda)));

        V_c += hbarc*(((1./(4.*pi*r))*coeffc[uPar]*(pow(pow(lambda,2) /
         (pow(lambda,2) - pow(mn[uPar],2)),2)))*(exp(-mn[uPar]*r)-exp(-lambda*r)*
        ((pow(lambda,2) - pow(mn[uPar],2))*r + 2*lambda) / (2*lambda)));
}
        Result =  V_rest + beta1 * V_c;
        return Result;
} 

//void AccessTObjArray() {
//    // Open the ROOT file
//    TFile* file = TFile::Open("Output_pOmega_Gamma_itResults.root", "READ");
//    if (!file || file->IsZombie()) {
//        std::cerr << "Error: Could not open the file." << std::endl;
//        return;
//    }
//
//    // List objects in the file for debugging
//    file->ls();
//  TCanvas *c2 = new TCanvas("c2", "p-#Omega with Betas", 800, 600);
//    // Loop through each TGraph by checking the keys in the file
//    for (int i = 1; i <= 100; ++i) {  // Assuming you want to check up to GammaGraphs;100
//        TString graphName = Form("GammaGraphs;%d", i);
//        TGraph* graph1 = dynamic_cast<TGraph*>(file->Get(graphName));
//
//        if (graph1) {
//            std::cout << "Retrieved TGraph " << graphName << " with " << graph1->GetN() << " points." << std::endl;
//            graph1->SetName(graphName); // Rename for identification if needed
//            //graph1->Draw("AL"); // Draw the graph on the current canvas
//        } else {
//            std::cerr << "Error: Could not retrieve " << graphName << "!" << std::endl;
//        }
//    }
//    c2->Update();
//    c2->SaveAs("betas.pdf");
//
//
//    // Close the file
//    file->Close();
//}

void readRootFile() {
// THE ERRORS ARE WRONG!!!! JUST A TEST TO SEE HOW TO PLOT THEM.
// Moving output data from remote ssh connection to local
// scp marta@elevation.lnf.infn.it:~/CATS_Xpi/Output_pOmega_Lattice.root ~/Desktop/Marta/Uni/MASTER/interns/INFN

    // Open the .root file
    TFile *datafile = TFile::Open("HEPData-ins1797617-v2-Table_2.root");
    
    if (!datafile || datafile->IsZombie()) {
        std::cerr << "Error opening file: " << "HEPData-ins1797617-v2-Table_2.root" << std::endl;
        return;
    }
    datafile->ls();
    
    // Access the directory "Table 2"
    TDirectoryFile *tableDir = (TDirectoryFile*)datafile->Get("Table 2");
    if (!tableDir) {
        std::cerr << "Directory 'Table 2' not found!" << std::endl;
        return;
    }
    
    // Access the original 2D histogram
    TH2F *histogram = (TH2F*)tableDir->Get("Hist2D_y1");
    if (!histogram) {
        std::cerr << "Histogram 'Hist2D_y1' not found!" << std::endl;
        return;
    }

     // Create arrays for the TGraph (X values and Z values)
    int nBinsX = histogram->GetNbinsX();
    int nBinsY = histogram->GetNbinsY();
    
    // Arrays to store X and Z values
    std::vector<double> xVals;
    std::vector<double> zVals;
    std::vector<double> zErrors;
    
    // Loop over the original histogram and fill X and Z values (Y -> Z)
    for (int i = 1; i <= nBinsX; i++) {
        for (int j = 1; j <= nBinsY; j++) {
            // Get the X and Y bin centers
            double xValue = histogram->GetXaxis()->GetBinCenter(i);
            double zValue = histogram->GetBinContent(i, j);  // Bin content (Z) becomes Y axis

             // Get the error associated with the Z value directly from the histogram
            double zError = histogram->GetBinError(i, j);  // Directly get the error from the histogram
            zError = zError/10;
            if (zValue != 0) {  // Only store points where Z has data
                xVals.push_back(xValue);
                zVals.push_back(zValue);
                zErrors.push_back(zError);
            }
        }
    }
 std::cout << "Number of data points: " << xVals.size() << std::endl;

    // Create a TGraph using X and Z values
   TGraphErrors *graph = new TGraphErrors(xVals.size(), &xVals[0], &zVals[0], nullptr, &zErrors[0]);
    graph->SetTitle("p-#Omega with Lattice Model; k^{*} (Mev/c); C(k^{*})");
    graph->SetMarkerStyle(kDot);
    graph->SetMarkerSize(4);





    // Open the ROOT file with all the betas
    TFile* file = TFile::Open("Output_pOmega_Gamma_itResults.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return;
    }

    // List objects in the file for debugging
    file->ls();


    //TCanvas *c2 = new TCanvas("c2", "p-#Omega with Betas", 800, 600);
    
    TGraph* graph1 = nullptr;
    
    // Loop through each TGraph by checking the keys in the file
    for (int i = 1; i <= 100; ++i) {  // Assuming you want to check up to GammaGraphs;100
        TString graphName = Form("GammaGraphs;%d", i);
        TGraph* tmpGraph = dynamic_cast<TGraph*>(file->Get(graphName));

        if (tmpGraph) {
            std::cout << "Retrieved TGraph " << graphName << " with " << tmpGraph->GetN() << " points." << std::endl;
            tmpGraph->SetName(graphName); // Rename for identification if needed
            
            // Check if this is the specific graph we want
            if (graphName == "GammaGraphs;3") {
                graph1 = tmpGraph; // Store it for later use
            }
        } else {
            std::cerr << "Error: Could not retrieve " << graphName << "!" << std::endl;
        }
    }

    //Open the second file (contains the model)

    TFile *modelFile = TFile::Open("Output_pOmega.root");
    if (!modelFile || modelFile->IsZombie()) {
        std::cerr << "Error opening model file: " << "Output_pOmega.root" << std::endl;
        return;
    }
    modelFile->ls();
    TGraph *modelGraph = (TGraph*)modelFile->Get("FitResult_pOmega");
    if (!modelGraph) {
        std::cerr << "Model graph 'FitResult_pOmega' not found!" << std::endl;
        return;
    }
    

    // Open the third file 
    TFile *modelFile2 = TFile::Open("Output_pOmega_Lattice_oton.root");
    if (!modelFile2 || modelFile2->IsZombie()) {
        std::cerr << "Error opening model file: " << "Output_pOmega_Lattice_oton.root" << std::endl;
        return;
    }
    
    modelFile2->ls();
    // Assume the model is stored as a TGraph or TH1, for example, "modelGraph"
    TGraph *modelGraph2 = (TGraph*)modelFile2->Get("FitResult_pOmega");
    if (!modelGraph2) {
        std::cerr << "Model graph2 'FitResult_pOmega' not found!" << std::endl;
        return;
    }

    // Open the fourth file 
    TFile *modelFile3 = TFile::Open("Output_pOmega_meson.root");
    if (!modelFile3 || modelFile3->IsZombie()) {
        std::cerr << "Error opening model file: " << "Output_pOmega_meson.root" << std::endl;
        return;
    }
    modelFile3->ls();
    // Assume the model is stored as a TGraph or TH1, for example, "modelGraph"
    TGraph *modelGraph3 = (TGraph*)modelFile3->Get("FitResult_pOmega");
    if (!modelGraph3) {
        std::cerr << "Model graph3 'FitResult_pOmega' not found!" << std::endl;
        return;
    }



    // Open the fifth file 
    TFile *modelFile4 = TFile::Open("Output_pOmega_meson_marta.root");
    if (!modelFile4 || modelFile4->IsZombie()) {
        std::cerr << "Error opening model file: " << "Output_pOmega_meson_marta.root" << std::endl;
        return;
    }
    modelFile4->ls();
    // Assume the model is stored as a TGraph or TH1, for example, "modelGraph"
    TGraph *modelGraph4 = (TGraph*)modelFile4->Get("FitResult_pOmega");
    if (!modelGraph4) {
        std::cerr << "Model graph4 'FitResult_pOmega' not found!" << std::endl;
        return;
    }

    // Open the sixt file 
    TFile *modelFile5 = TFile::Open("Output_pOmega_meson_0.3000000000.root");
    if (!modelFile5 || modelFile5->IsZombie()) {
        std::cerr << "Error opening model file: " << "Output_pOmega_meson_0.3000000000.root" << std::endl;
        return;
    }
    modelFile5->ls();
    // Assume the model is stored as a TGraph or TH1, for example, "modelGraph"
    TGraph *modelGraph5 = (TGraph*)modelFile5->Get("FitResult_pOmega");
    if (!modelGraph5) {
        std::cerr << "Model graph4 'FitResult_pOmega' not found!" << std::endl;
        return;
    }
    // Open the seventh file 
    TFile *modelFile6 = TFile::Open("Output_pOmega_meson_0.03.root");
    if (!modelFile6 || modelFile6->IsZombie()) {
        std::cerr << "Error opening model file: " << "Output_pOmega_meson_0.03.root" << std::endl;
        return;
    }
    modelFile6->ls();
    // Assume the model is stored as a TGraph or TH1, for example, "modelGraph"
    TGraph *modelGraph6 = (TGraph*)modelFile6->Get("FitResult_pOmega");
    if (!modelGraph6) {
        std::cerr << "Model graph4 'FitResult_pOmega' not found!" << std::endl;
        return;
    }
        // Open the seventh file 
    TFile *modelFile7 = TFile::Open("Output_pOmega_meson_0.1000000000.root");
    if (!modelFile7 || modelFile7->IsZombie()) {
        std::cerr << "Error opening model file: " << "Output_pOmega_meson_0.1000000000.root" << std::endl;
        return;
    }
    modelFile7->ls();
    // Assume the model is stored as a TGraph or TH1, for example, "modelGraph"
    TGraph *modelGraph7 = (TGraph*)modelFile7->Get("FitResult_pOmega");
    if (!modelGraph7) {
        std::cerr << "Model graph4 'FitResult_pOmega' not found!" << std::endl;
        return;
    }

/////// FIT FUNCTION
TF1 *Tuning = new TF1("meson tuning", Tetsuo_tuning,-100, 400, 1); // Example custom function
//Tuning->SetParameters(1, 1, 1); // Initial parameter guesses



    // Step 3: Create a canvas and plot the data first
    TCanvas *c1 = new TCanvas("c1", "p-#Omega with Lattice Model", 800, 600);
    graph->Draw("AP");

    // Perform the fit and draw it on the same canvas
    //graph->Fit(Tuning, "R"); // Fit the data with the custom function
    //Tuning->Print();
    //std::cout << "Chi2: " << Tuning->GetChisquare() << std::endl;
    //std::cout << "Ndf: " << Tuning->GetNDF() << std::endl; // This should give you the degrees 

    // Step 4: Plot the model on top of the data
    modelGraph->SetLineColor(kGreen);  // Set model color to red
    modelGraph->SetLineWidth(3);     // Set model line width
    modelGraph->Draw("L SAME");      // Draw model as a line (L) on the same canvas

    // Plot the second model on top of the data
    modelGraph2->SetLineColor(kBlue);  // Set model color to red
    modelGraph2->SetLineStyle(4);
    modelGraph2->SetLineWidth(3);     // Set model line width
    modelGraph2->Draw("L SAME");      // Draw model as a line (L) on the same canvas

    // Plot the third model on top of the data
    modelGraph3->SetLineColor(9);  // Set model color to red
    modelGraph3->SetLineStyle(1);
    modelGraph3->SetLineWidth(3);     // Set model line width
    modelGraph3->Draw("L SAME");      // Draw model as a line (L) on the same canvas


    // Plot the fourth model on top of the data
    modelGraph4->SetLineColor(7);  // Set model color to cyan
    modelGraph4->SetLineStyle(4);
    modelGraph4->SetLineWidth(3);     // Set model line width
    modelGraph4->Draw("L SAME");      // Draw model as a line (L) on the same canvas

    // Plot the fourth model on top of the data
    modelGraph5->SetLineColor(9);  // Set model color to cyan
    modelGraph5->SetLineStyle(1);
    modelGraph5->SetLineWidth(3);     // Set model line width
    modelGraph5->Draw("L SAME");      // Draw model as a line (L) on the same canvas
    // Plot the fourth model on top of the data
    //modelGraph6->SetLineColor(28);  // Set model color to purple
    //modelGraph6->SetLineStyle(1);
    //modelGraph6->SetLineWidth(3);     // Set model line width
    //modelGraph6->Draw("L SAME");      // Draw model as a line (L) on the same canvas
    // Plot the fourth model on top of the data
    //modelGraph7->SetLineColor(2);  // Set model color to purple
    //modelGraph7->SetLineStyle(1);
    //modelGraph7->SetLineWidth(3);     // Set model line width
    //modelGraph7->Draw("L SAME");      // Draw model as a line (L) on the same canvas


    // Check if the graph was found
    if (graph1) {
        graph1->SetMarkerSize(4);
        graph1->SetLineColor(kRed); // Optionally set color for better visibility
        graph1->Draw("L SAME"); // Draw the graph on the canvas
    } else {
        std::cerr << "GammaGraphs;3 was not found!" << std::endl;
    }


    // Step 5: Add a legend
    TLegend *legend = new TLegend(0.4, 0.67, 0.9, 0.9);
    legend->AddEntry(graph, "ALICE data", "p");
    legend->AddEntry(graph1, " tuning", "l");
    legend->AddEntry(modelGraph, "Lattice QCD Model Marta", "l");
    legend->AddEntry(modelGraph2, "Lattice QCD Model Oton", "l");
    legend->AddEntry(modelGraph3, "meson exchange model CATS", "l");
    legend->AddEntry(modelGraph4, "meson exchange model Marta", "l");
    legend->AddEntry(modelGraph5, "meson exchange model #beta = 0.3", "l");
    //legend->AddEntry(modelGraph6, "meson exchange model #beta = 0.03", "l");
    //legend->AddEntry(modelGraph7, "meson exchange model #beta = 0.1", "l");
    //legend->AddEntry(Tuning, "Fit: Custom function", "l");  // Ad
   // legend->AddEntry(Tuning, "meson exchange model #beta = 0.0001", "l");
    legend->SetTextSize(0.03); 
    legend->Draw();

    // Step 6: Save and update the canvas
    c1->Update();
    c1->SaveAs("data_with_model.pdf");



    // Close the file
    file->Close();
   

    // Close files
    datafile->Close();
    modelFile->Close();   
    modelFile2->Close(); 
    modelFile3->Close(); 
    modelFile4->Close(); 
    modelFile5->Close(); 
    modelFile6->Close(); 
}


int main(){
    readRootFile();
    }