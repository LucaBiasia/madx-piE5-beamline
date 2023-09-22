// apro il file e prendo il TTree
void TreeToNTuple() {
    TFile* f = new TFile("/home/lucabiasia/Scrivania/single_element/Init_tmp.root");
    TTree* init = (TTree*)f->Get("Init_tmp");

    TFile *file = new TFile("InitPhaseSpace.root", "RECREATE");
    TNtuple *starter = new TNtuple("ntuple", "inizialize_ps", "x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight");
    if(!init){
        // il file non contiene il TTree init
        printf("cannot find the tree\n");
        return;
    }
    // init->Scan();
    float x, y, z, Px, Py, Pz, t, PDGid, EventID, TrackID, ParentID, Weight;
    init->SetBranchAddress("x", &x);
    init->SetBranchAddress("y", &y);
    init->SetBranchAddress("z", &z);
    init->SetBranchAddress("Px", &Px);
    init->SetBranchAddress("Py", &Py);
    init->SetBranchAddress("Pz", &Pz);
    init->SetBranchAddress("t", &t);
    init->SetBranchAddress("PDGid", &PDGid);
    init->SetBranchAddress("EventID", &EventID);
    init->SetBranchAddress("TrackID", &TrackID);
    init->SetBranchAddress("ParentID", &ParentID);
    init->SetBranchAddress("Weight", &Weight);
    
    Int_t nEntries = init->GetEntries();
    // printf("%d eventi\n", nEntries);

    for(Int_t iEntry=0; iEntry < nEntries; iEntry++){
        init->GetEntry(iEntry);
        starter->Fill(x, y, z, Px, Py, Pz, t, PDGid, EventID, TrackID, ParentID, Weight);

    }
    file->cd();
    starter->Write();

    
    // starter->Draw("x:y");
    f->Close();
    file->Close();
}
