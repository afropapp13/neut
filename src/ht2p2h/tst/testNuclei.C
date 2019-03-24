{ 
  Nucleus2p2h nn;

  nn.Initialize(12);
  
  TH1F *hfm = new TH1F("hfm","   ",100,0.,7.0);

  TH1F *hr = new TH1F("hr","   ",100,0.,7.0); 

  TH1F *hr2 = new TH1F("hr2","   ",100,0.,7.0); 

  double tot = 0.; 
  
  for(int i = 0 ; i < 100; i++ ) {
    
    hr->SetBinContent(i+1,nn.Density(12,1,hr->GetBinCenter(i+1)));

    tot += nn.Density(12,1,hr->GetBinCenter(i+1));
    
    hfm->SetBinContent(i+1,nn.GetFermiLFG(hr->GetBinCenter(i+1),12,1)); 
  }

  for( int i = 0; i < 100000; i++  ) {    
    hr2->Fill(nn.GenerateR(12,1),1./100000.*tot); 
  }
  
  hfm->Draw();
  hr->Draw("same");
  hr2->Draw("same"); 
  
}
