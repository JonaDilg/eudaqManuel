#include "CaribouEvent2StdEventConverter.hh"
#include "utils/log.hpp"
#include "CLICTDFrameDecoder.hpp"
#include "TFile.h"
#include "TH1D.h"
#include "time.h"

using namespace eudaq;

namespace{
    auto dummy0 = eudaq::Factory<eudaq::StdEventConverter>::
        Register<DSO9254AEvent2StdEventConverter>(DSO9254AEvent2StdEventConverter::m_id_factory);
}

// i use these
bool    DSO9254AEvent2StdEventConverter::m_configured(0);

bool    DSO9254AEvent2StdEventConverter::m_generateRoot(0);
bool DSO9254AEvent2StdEventConverter::m_osci_timestamp(1);


// not used in reduced version, see .bak file for full version
int64_t DSO9254AEvent2StdEventConverter::m_runStartTime(-1);
double DSO9254AEvent2StdEventConverter::m_pedStartTime(0);
double DSO9254AEvent2StdEventConverter::m_pedEndTime(0);
double DSO9254AEvent2StdEventConverter::m_ampStartTime(0);
double DSO9254AEvent2StdEventConverter::m_ampEndTime(0);
double DSO9254AEvent2StdEventConverter::m_chargeScale(0);
double DSO9254AEvent2StdEventConverter::m_chargeCut(0);

bool DSO9254AEvent2StdEventConverter::m_polarity(1);

bool DSO9254AEvent2StdEventConverter::Converting(eudaq::EventSPC d1, eudaq::StandardEventSP d2, eudaq::ConfigurationSPC conf) const{

    auto event = std::dynamic_pointer_cast<const eudaq::RawEvent>(d1);
    // No event
    if(!event)
        return false;
    
    std::ofstream outfileTimestamps;
    if( !m_configured ) {
        // read from config file
        m_generateRoot = conf->Get("generateRoot", 0);
        m_osci_timestamp = conf->Get("osci_timestamp", 1);

        // check configuration
        bool noData = false;
        bool noHist = false;
        if( !m_generateRoot ){
            EUDAQ_WARN( "No waveform file generated!" );
            noHist = true;
        }

        m_configured = true;
    } // configure

    // Data container:
    caribou::pearyRawData rawdata;

    // internal containers
    std::vector<std::vector<std::vector<double>>> wavess;
    std::vector<double> time;
    uint64_t timestamp;

    // waveform parameters
    std::vector<uint64_t> np;
    std::vector<double> dx;
    std::vector<double> x0;
    std::vector<double> dy;
    std::vector<double> y0;

    // Bad event
    if(event->NumBlocks() != 1) {
        EUDAQ_WARN("Ignoring bad packet " + std::to_string(event->GetEventNumber()));
        return false;
    }

    // generate rootfile to write waveforms as TH1D
    TFile * histoFile = nullptr;
    if( m_generateRoot ){
        histoFile = new TFile( Form( "waveforms_run%i.root", event->GetRunN() ), "UPDATE" );
        if(!histoFile) {
            EUDAQ_ERROR(to_string(histoFile->GetName()) + " can not be opened");
            return false;
        }
    }


    /* all four scope channels in one data block */
    auto datablock = event->GetBlock(0);
    EUDAQ_DEBUG("DSO9254A frame with " + to_string(datablock.size()) + " entries");

    // Calulate positions and length of data blocks:
    // FIXME FIXME FIXME by Simon: this is prone to break since you are selecting bits from a 64bit
    //                             word - but there is no guarantee in the endianness here and you
    //                             might end up having the wrong one.
    // Finn: I am not sure if I actually solved this issue! Lets discuss.
    rawdata.resize( sizeof(datablock[0]) * datablock.size() / sizeof(uintptr_t) );
    std::memcpy(&rawdata[0], &datablock[0], sizeof(datablock[0]) * datablock.size() );

    // needed per event
    uint64_t block_words;
    uint64_t pream_words;
    uint64_t chann_words;
    uint64_t block_position = 0;
    uint64_t sgmnt_count = 1;

    // loop 4 channels
    for( int i_channel = 0; i_channel < 4; i_channel++ ){
        EUDAQ_DEBUG( "Reading channel " + to_string(i_channel));

        // Obtaining Data Stucture
        block_words = rawdata[block_position + 0];
        pream_words = rawdata[block_position + 1];
        chann_words = rawdata[block_position + 2 + pream_words];
        EUDAQ_DEBUG( "  " + to_string(block_position) + " current position in block");
        EUDAQ_DEBUG( "  " + to_string(block_words) + " words per segment");
        EUDAQ_DEBUG( "  " + to_string(pream_words) + " words per preamble");
        EUDAQ_DEBUG( "  " + to_string(chann_words) + " words per channel data");
        EUDAQ_DEBUG( "  " + to_string(rawdata.size()) + " size of rawdata");

        // Check our own sixth-graders math: The block header minus peamble minus channel data minus the counters is zero
        if(block_words - pream_words - 2 - chann_words) {
            EUDAQ_WARN("Go back to school, 6th grade - math doesn't check out! :/");
        }

        // read preamble:
        std::string preamble = "";
        for(int i=block_position + 2;
                i<block_position + 2 + pream_words;
                i++ ){
            preamble += (char)rawdata[i];
        }
        EUDAQ_DEBUG("Preamble");
        EUDAQ_DEBUG("  " + preamble);
        // parse preamble to vector of strings
        std::string val;
        std::vector<std::string> vals;
        std::istringstream stream(preamble);
        while(std::getline(stream,val,',')){
            vals.push_back(val);
            EUDAQ_DEBUG( " " + to_string(vals.size()-1) + ": " + to_string(vals.back()));
        }

        // Pick needed preamble elements
        np.push_back( stoi( vals[2]) );
        dx.push_back( stod( vals[4]) * 1e9 );
        x0.push_back( stod( vals[5]) * 1e9 );
        dy.push_back( stod( vals[7]) );
        y0.push_back( stod( vals[8]) );

        EUDAQ_DEBUG("Acq mode (2=SEGM): " + to_string(vals[18]));
        if( vals.size() == 25 ) {// this is segmented mode, possibly more than one waveform in block
            EUDAQ_DEBUG( "Segments: " + to_string(vals[24]));
            sgmnt_count = stoi( vals[24] );
        }

        // Check our own seventh-graders math: total data size minus four times the data of a channel minus four block header words = 0
        if(rawdata.size() - 4 * block_words - 4) {
            EUDAQ_WARN("Go back to school 7th grade - math doesn't check out! :/ " + to_string(rawdata.size() - 4 * block_words) + " is not 0!");
        }

        int points_per_words = np.at(i_channel)/(chann_words/sgmnt_count);
        if(chann_words % sgmnt_count != 0) {
            EUDAQ_WARN("Segment count and channel words don't match: " + to_string(chann_words) + "/" + to_string(sgmnt_count));
        }
        if( np.at(i_channel) % (chann_words/sgmnt_count) ){ // check
            EUDAQ_WARN("incomplete waveform in block " + to_string(event->GetEventN())
                                 + ", channel " + to_string(i_channel) );
        }
        EUDAQ_DEBUG("WF points per data word: " + to_string(points_per_words));

        // Check our own eighth-graders math: the number of words in the channel times 4 points per word minus the number of points times number of segments is 0
        if(4 * chann_words - np.at(i_channel) * sgmnt_count) {
            EUDAQ_WARN("Go back to school 8th grade - math doesn't check out! :/ " + to_string(chann_words - np.at(i_channel) * sgmnt_count) + " is not 0!");
        }

        // set run start time to 0 ms, using  timestamps from preamble which are 10 ms precision
        if( m_runStartTime < 0 ){
            m_runStartTime = DSO9254AEvent2StdEventConverter::timeConverter( vals[15], vals[16] );
        }
        timestamp = ( DSO9254AEvent2StdEventConverter::timeConverter( vals[15], vals[16] )
                                    - m_runStartTime ) * 1e9; // ms to ps

        // not use osc timestamp:
        if( !m_osci_timestamp ){
            timestamp = 0;
        }

        // need once per channel
        std::vector<double> ped;
        std::vector<double> amp;
        std::vector<std::vector<double>> waves;

        // loop over segments
        for (int i_segment = 0; i_segment < sgmnt_count; i_segment++) {
            // contains one waveform:
            std::vector<double> wave;

            TH1D * hist = nullptr;
            if( m_generateRoot ){
                TH1D* hist_init = new TH1D( Form( "waveform_run%i_ev%i_ch%i_s%i", event->GetRunN(), event->GetEventN(), i_channel, i_segment ),
                                            Form( "waveform_run%i_ev%i_ch%i_s%i", event->GetRunN(), event->GetEventN(), i_channel, i_segment ),
                                            np.at(i_channel),
                                            x0.at(i_channel) - dx.at(i_channel)/2.,
                                            x0.at(i_channel) - dx.at(i_channel)/2. + dx.at(i_channel)*np.at(i_channel));
                hist = hist_init;
                hist->GetXaxis()->SetTitle("time [ns]");
                hist->GetYaxis()->SetTitle("signal [V]");
            }

            // read channel data
            std::vector<int16_t> words;
            words.resize(points_per_words); // rawdata contains 8 byte words, scope sends 2 byte words
            int16_t i_waveform = 0;
            // Read from segment start until the next segment begins:
            for(int i=block_position + 3 + pream_words + (i_segment+0)*chann_words/sgmnt_count;
                    i<block_position + 3 + pream_words + (i_segment+1)*chann_words/sgmnt_count;
                    i++){

                // copy channel data from entire segment data block
                memcpy(&words.front(), &rawdata[i], points_per_words*sizeof(int16_t));

                for( auto word : words ){

                    // fill vectors with time bins and waveform
                    if( i_channel == 0 && i_segment == 0 ){ // this we need only once
                        time.push_back( i_waveform * dx.at(i_channel) + x0.at(i_channel) );
                    }
                    wave.push_back( (word * dy.at(i_channel) + y0.at(i_channel)) );


                    // fill histogram
                    if( m_generateRoot ){
                        hist->SetBinContent( hist->FindBin( i_waveform * dx.at(i_channel) + x0.at(i_channel) ),
                                                                 (double)word * dy.at(i_channel) + y0.at(i_channel) );
                    }

                    i_waveform++;

                } // for: words

            } // for: waveform

            // store waveform vector and create vector entries for pedestal and amplitude
            waves.push_back( wave );

            // write and delete
            if( m_generateRoot ){
                hist->Write();
            }
            delete hist;

        } // for: i_segment

        wavess.push_back( waves );

        // update position for next iteration
        block_position += block_words+1;
    } // for: i_channel

    // close rootfile
    if( m_generateRoot ){
        histoFile->Close();
        EUDAQ_DEBUG("Histograms written to " + to_string(histoFile->GetName()));
    }

    // declare map between scope channel number and 2D pixel index
    // TODO make this configurable
    std::map<int,std::vector<int>> chanToPix;
    chanToPix[0] = {1,0};
    chanToPix[1] = {0,0};
    chanToPix[2] = {1,1};
    chanToPix[3] = {0,1};

    // process waveforms
    if(time.size() > 0 ){      // only if there are waveform data

        // re-iterate waveforms
        // wavess.at(channel).at(segment).at(point)
        // ped   .at(channel).at(segment)
        for( int i_segment = 0; i_segment<wavess.at(0).size(); i_segment++ ){
            for( int i_channel = 0; i_channel<wavess.size(); i_channel++ ){
                /* Jona: in old implementation, pedestal and amplitude would be calcualted here*/

                // Checking:
                EUDAQ_DEBUG("CHECK: channel "+ to_string(i_channel) + " segment " + to_string(i_segment));
                EUDAQ_DEBUG("  points  " + to_string(wavess.at(i_channel).at(i_segment).size()));
                EUDAQ_DEBUG("  timestamp " + to_string(timestamp));
            } // for: i_channel

            // create sub-event for segments
            auto sub_event = eudaq::StandardEvent::MakeShared();

            // Create a StandardPlane representing one sensor plane
            eudaq::StandardPlane plane(0, "Caribou", "DSO9254A");

             // index number, which is the number of pixels hit per event
            int index = 0;

            // fill plane
            plane.SetSizeZS(2, 2, 0);
            for( int i_channel = 0; i_channel<wavess.size(); i_channel++ ){
                /* Jona: here, only pixels above threshold were selected. i take them all. */
                plane.PushPixel( chanToPix[i_channel].at(0), chanToPix[i_channel].at(1), 0, timestamp );

                // Set waveforms to each hit pixel.
                plane.SetWaveform(index, wavess.at(i_channel).at(i_segment), time.at(i_channel), dx.at(i_channel), 0);

                index++;
            } // for: i_channel

            // derive trigger number from block number
            int triggerN = (event->GetEventN()-1) * wavess.at(0).size() + i_segment + 1;

            EUDAQ_DEBUG("Block number " + to_string(event->GetEventN()) + " " +
                                    " block size " + to_string(wavess.at(0).size()) +
                                    " segments, segment number " + to_string(i_segment) +
                                    " trigger number " + to_string(triggerN));

            sub_event->AddPlane(plane);
            sub_event->SetDetectorType("DSO9254A");
            sub_event->SetTimeBegin(0); // forcing corry to fall back on trigger IDs
            sub_event->SetTimeEnd(0);
            sub_event->SetTriggerN(triggerN);

            // add sub event to StandardEventSP for output to corry
            d2->AddSubEvent( sub_event );

        } // for: i_segment

        EUDAQ_DEBUG( "Number of sub events " + to_string(d2->GetNumSubEvent()) );

    } // if there is waveform data
    else {
        EUDAQ_WARN("No scope data in block " + to_string(event->GetEventN()));
        return false;
    }

    // Indicate that data were successfully converted
    return true;
}



uint64_t DSO9254AEvent2StdEventConverter::timeConverter( std::string date, std::string time ){

    // needed for month conversion
    std::map<std::string, int> monthToNum {
        {"JAN",  1},
        {"FEB",  2},
        {"MAR",  3},
        {"APR",  4},
        {"MAI",  5},
        {"JUN",  6},
        {"JUL",  7},
        {"AUG",  8},
        {"SEP",  9},
        {"OCT", 10},
        {"NOV", 11},
        {"DEC", 12},
    };

    // delete leading and tailing "
    date.erase( 0, 1);
    date.erase( date.size()-1, 1 );
    time.erase( 0, 1);
    time.erase( time.size()-1, 1 );

    // get day, month ... in one vector of strings
    std::string s;
    std::istringstream s_date( date );
    std::istringstream s_time( time );
    std::vector<std::string> parts;
    while(std::getline(s_date,s,' ')) parts.push_back(s);
    while(std::getline(s_time,s,':')) parts.push_back(s);

    // fill time structure
    struct std::tm build_time = {0};
    build_time.tm_mday = stoi( parts.at(0) );       // day in month 1-31
    build_time.tm_year = stoi( parts.at(2) )-1900;  // years since 1900
    build_time.tm_hour = stoi( parts.at(3) );       // hours since midnight 0-23
    build_time.tm_min  = stoi( parts.at(4) );       // minutes 0-59
    build_time.tm_sec  = stoi( parts.at(5) );       // seconds 0-59
    build_time.tm_mon  = monthToNum[parts.at(1)]-1; // month 0-11
    if( monthToNum[parts.at(1)]==0 ){ // failed conversion
        EUDAQ_ERROR(" in timeConverter, month " + to_string(parts.at(1)) + " not defined!");
    }

    // now convert to unix timestamp, to milli seconds and add sub-second information from scope
    std::time_t tunix;
    tunix = std::mktime( &build_time );
    uint64_t result = static_cast<uint64_t>( tunix ) * 1e3;               // [   s->ms]
    uint64_t mill10 = static_cast<uint64_t>( stoi( parts.at(6) ) ) * 1e1; // [10ms->ms]
    result += mill10;

    return result;
} // timeConverter
