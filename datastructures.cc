// Datastructures.cc
//
// Student name: Mikko Kaukonen


#include "datastructures.hh"

#include <random>

#include <cmath>

#include <map>
#include <vector>
#include <set>
#include <unordered_set>
#include <iterator>
#include <algorithm>
#include <string>
#include <utility>
#include <iostream>
#include <stack>
#include <list>


// This function was ready in the program code base of the exercise.
// So not done by me!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

std::minstd_rand rand_engine; // Reasonably quick pseudo-random generator

template <typename Type>
Type random_in_range(Type start, Type end)
{
    auto range = end-start;
    ++range;

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine);

    return static_cast<Type>(start+num);
}

// Modify the code below to implement the functionality of the class.
// Also remove comments from the parameter names when you implement
// an operation (Commenting out parameter name prevents compiler from
// warning about unused parameters on operations you haven't yet implemented.)

Datastructures::Datastructures()
{
    // Write any initialization you need here
}

/*!
 * \brief Datastructures::~Datastructures tyhjentää tietorakenteet
 */
Datastructures::~Datastructures()
{
    clear_all();

}

/*!
 * \brief Datastructures::station_count Laskee asemien määrän tietorakenteesta
 * \return asemien lukumäärä
 */
unsigned int Datastructures::station_count()
{
    return stationdata_.size();
}

/*!
 * \brief Datastructures::clear_all Tyhjentää kaikki luokan hh private osioon
 * luodut tietorakenteet.
 */
void Datastructures::clear_all()
{
    {
        {
        for(std::unordered_map<std::string, Station*>::iterator station = stationdata_.begin();
            station != stationdata_.end(); station++)
            {
            delete (station->second);
            }
        stationdata_.clear();
        }

        {
        for(std::unordered_map<unsigned long long int, Region*>::iterator region = regionsdata_.begin();
            region != regionsdata_.end(); region++)
            {
            delete (region->second);
            }
        regionsdata_.clear();
        }


        sortingstation_.clear();

        clear_trains();

    }

}

/*!
 * \brief Datastructures::all_stations Palauttaa kaikki tietorakenteessa
 * olevat asemat mielivaltaisessa järjestyksessä
 * \return vektorin joka sisältää asemien id numerot
 */
std::vector<StationID> Datastructures::all_stations()
{
    std::vector<std::string> station_names;
    for(auto & entry :stationdata_)
        station_names.push_back(entry.first);
    return station_names;

}

/*!
 * \brief Datastructures::calculate_distance Laskee kahden koordinaatin
 * etäisyyden toisistaan.
 * \param x x koordinaatti
 * \param y y koordinaatti
 * \return koordinattien etäisyys toisistaan
 */
int Datastructures::calculate_distance(int x, int y)
{
    int Distance = sqrt(pow(x, 2) + pow(y, 2));
    return Distance;
}

/*!
 * \brief Datastructures::add_station Lisää tietorakenteeseen uuden aseman
 * annetulla uniikilla id:llä, nimellä ja sijainnilla.
 * \param id aseman id numero
 * \param name aseman nimi
 * \param xy aseman koordinaatit
 * \return totuusarvon true tai false
 */
bool Datastructures::add_station(StationID id, const Name& name, Coord xy)
{
    std::map<TrainID, Train*> trains;

    if(stationdata_.find(id) == stationdata_.end()){

        std::map<unsigned short int,std::vector<std::string>> depatures;

        Station* new_station = new Station{id, name, xy,
                calculate_distance(xy.x, xy.y), xy.y, xy.x, depatures,
                NO_REGION, trains, -1, nullptr, "white", {NO_STATION}, 0 };

        stationdata_.insert({id, new_station});
        sortingstation_.insert({name, new_station});
        return true;

    }

    else
        return false;
}

/*!
 * \brief Datastructures::get_station_name Palauttaa annetulla ID:llä olevan
 * aseman nimen.
 * \param id aseman id numero
 * \return aseman nimi
 */
Name Datastructures::get_station_name(StationID id)
{
    if(stationdata_.find(id) != stationdata_.end()){
        return stationdata_[id]->name;
    }
    else
        return NO_NAME;
}

/*!
 * \brief Datastructures::get_station_coordinates Palauttaa annetulla ID:llä
 * olevan aseman koordinaatit.
 * \param id aseman id numero
 * \return aseman koordinaatit
 */
Coord Datastructures::get_station_coordinates(StationID id)
{
    if(stationdata_.find(id) != stationdata_.end()){
        return stationdata_[id]->coord;
    }
    else
        return NO_COORD;
}

/*!
 * \brief Datastructures::stations_alphabetically Palauttaa id:t asemien
 * nimen mukaan aakkosjärjestyksessä.
 * \return Vektori jossa asemien id:t järjestyksessä
 */
std::vector<StationID> Datastructures::stations_alphabetically()
{
    std::vector<std::string> sortedstations;
    for(auto & entry : sortingstation_){
        sortedstations.push_back(entry.second->stationid);
    }

    return sortedstations;

}

/*!
 * \brief Datastructures::stations_distance_increasing Palauttaa id:t
 * asemien koordinaattien mukaan järjestettynä
 * \return Vektori jossa asemien id:t järjestyksessä
 */
std::vector<StationID> Datastructures::stations_distance_increasing()
{
    std::vector<Station*> vec;
    std::vector<std::string> sortedcoords;


    for(auto & entry : stationdata_){
        vec.push_back(entry.second);
    }
    std::sort(vec.begin(), vec.end(), [] (Station* a, Station* b)
    {
        if(a->Distance != b->Distance)
            {
            return a->Distance < b->Distance;
            }
        else {
            return a->y_coord < b->y_coord;}
            }
    );

    for(std::vector<Station*>::iterator it = vec.begin(); it != vec.end(); ++it)
        sortedcoords.push_back((*it)->stationid);
    return sortedcoords;
}

/*!
 * \brief Datastructures::find_station_with_coord Palauttaa aseman, joka on
 * annetussa koordinaatissa tai NO_STATION, jos sellaista ei ole.
 * \param xy struct joka sisältää x ja y koordinaatit
 * \return annetuissa koordinaateissa olevan aseman id numero
 */
StationID Datastructures::find_station_with_coord(Coord xy)
{
    for(auto & entry : stationdata_){
        if (entry.second->coord == xy)
            return entry.second->stationid;
    }
    return NO_STATION;
}

/*!
 * \brief Datastructures::change_station_coord Muuttaa annetulla ID:llä
 * olevan aseman sijainnin.
 * \param id aseman id numero
 * \param newcoord uudet koordinaatit
 * \return totuusarvo true tai false
 */
bool Datastructures::change_station_coord(StationID id, Coord newcoord)
{
    if(stationdata_.find(id) != stationdata_.end()){
        stationdata_[id]->coord=newcoord;
        stationdata_[id]->Distance=(calculate_distance(newcoord.x, newcoord.y));
        stationdata_[id]->y_coord=(newcoord.y);
        return true;
        }
    else
        return false;
}

/*!
 * \brief Datastructures::add_departure Lisää tiedon siitä, että annetulta
 * asemalta lähtee annettu juna annetulla kellonajalla.
 * \param stationid aseman id numero
 * \param trainid junan id numero
 * \param time junan lähtöaika
 * \return totuusarvo true tai false
 */
bool Datastructures::add_departure(StationID stationid, TrainID trainid, Time time)
{
    std::vector<TrainID> train;

    if(stationdata_.find(stationid) != stationdata_.end()){

        Station* & A = stationdata_.at(stationid);

        if(A->depatures.find(time) == A->depatures.end())
       {
            train.push_back(trainid);
            A->depatures.insert({time,train});
            return true;
            }

        else if((std::find(A->depatures[time].begin(),
                            A->depatures[time].end(), trainid) ==
                            A->depatures[time].end())){
                            A->depatures[time].push_back(trainid);
            return true;
            }
        else
            return false;
            }
    else
        return false;
}

/*!
 * \brief Datastructures::remove_departure Poistaa annetulta asemalta annetun
 * junan lähdön annetulla kellonajalla.
 * \param stationid aseman id numero
 * \param trainid junan id numero
 * \param time junan lähtöaika
 * \return totuusarvo true tai false
 */
bool Datastructures::remove_departure(StationID stationid, TrainID trainid, Time time)
{

    if(stationdata_.find(stationid) != stationdata_.end())
    {

        Station* & A = stationdata_[stationid];

        auto remove_time = A->depatures.find(time);

        if((remove_time != A->depatures.end())
                    &(std::find(A->depatures[time].begin(),
                                A->depatures[time].end(), trainid) !=
                                A->depatures[time].end()))
        {
            std::vector<std::string>::iterator remove_train = std::find(A->
                                            depatures[time].begin(), A->
                                            depatures[time].end(), trainid);
            if(remove_train == A->depatures[time].end())
                return false;
            else
                A->depatures[time].erase(remove_train);
            return true;
        }

        return true;

    }
    return false;
}

/*!
 * \brief sortbydepatures Palauttaa asemat kellonajan ja sen jälkeen junan
 * id numeron mukaisessa järjestyksessä
 * \param a Station struct
 * \param b Station struct
 * \return  asemat annetussa järjestyksessä
 */
bool sortbydepatures(std::pair<Time, TrainID> &a,
              std::pair<Time, TrainID> &b)
{
    if(a.first != b.first)
        return (a.first < b.first);
    else
        return (a.second < b.second);
}

/*!
 * \brief Datastructures::station_departures_after Listaa kaikki lähdöt
 * annetulta asemalta annettunakellonaikana tai sen jälkeen. Lähdöt listataan
 * kellonaikajärjestyksessä ja samalla kellonlyömällä lähtevät juna-id:n mukaan
 * \param stationid aseman id numero
 * \param time junan lähtöaika
 * \return Vektorin pareja jossa junien lähdöt annetulta asemalta on listattuna
 * kuvauksessa määriteltyyn järjestykseen.
 */
std::vector<std::pair<Time, TrainID>> Datastructures::station_departures_after(StationID stationid, Time time)
{

    std::vector<std::pair<Time, TrainID>> vectoreturn;
    vectoreturn.clear();

    if(stationdata_.find(stationid) != stationdata_.end())
    {
        Station* & A = stationdata_[stationid];
        for(auto & entry : A->depatures)
        {
            if(time <= entry.first)
            {
                for(unsigned int i=0; i<entry.second.size(); i++)
                    vectoreturn.push_back(make_pair(entry.first, entry.second.at(i)));

            }

        }
        std::sort(vectoreturn.begin(), vectoreturn.end(), sortbydepatures);
        return vectoreturn;
    }
    vectoreturn.push_back({NO_TIME, NO_TRAIN});
    return vectoreturn;
}

/*!
 * \brief Datastructures::add_region Lisää tietorakenteeseen uuden alueen
 * annetulla uniikilla id:llä, nimellä ja monikulmiolla (koordinaatit)
 * \param id alueen id numero
 * \param name alueen nimi
 * \param coords alueen koordinaatit
 * \return totuusarvo true tai false
 */
bool Datastructures::add_region(RegionID id, const Name &name, std::vector<Coord> coords)
{
    if(regionsdata_.find(id) == regionsdata_.end())
    {
        std::vector<Coord> coordinates;
        std::vector<StationID> stations;
        std::set<Region*> children;

        Region* new_region = new Region{id, name, coords, stations, children, nullptr};
        regionsdata_.insert({id, new_region});
    return true;

    }
    else
        return false;
}

/*!
 * \brief Datastructures::all_regions Palauttaa kaikki tietorakenteessa olevat
 * alueet mielivaltaisessa järjestyksessä
 * \return Vektorin jossa alueiden id numerot
 */
std::vector<RegionID> Datastructures::all_regions()
{
    std::vector<RegionID> regions;
    for(auto & entry: regionsdata_)
        regions.push_back(entry.first);
    return regions;
}

/*!
 * \brief Datastructures::get_region_name Palauttaa annetulla ID:llä olevan
 * alueen nimen, tai NO_NAME, jos id:llä ei löydy aluetta.
 * \param id alueen id numero
 * \return alueen nimi
 */
Name Datastructures::get_region_name(RegionID id)
{
    if(regionsdata_.find(id) != regionsdata_.end())
        return regionsdata_.at(id)->name;
    else
        return (NO_NAME);
}

/*!
 * \brief Datastructures::get_region_coords Palauttaa annetulla ID:llä olevan
 * alueen koordinaattivektorin, tai vektorin, jonka ainoa alkio on NO_COORD,
 * jos id:llä ei löydy aluetta.
 * \param id alueen id numero
 * \return alueen koordinaatit
 */
std::vector<Coord> Datastructures::get_region_coords(RegionID id)
{
    if(regionsdata_.find(id) != regionsdata_.end())
        return regionsdata_[id]->coordinates;
    else
        return {NO_COORD};
}

/*!
 * \brief Datastructures::add_subregion_to_region Lisää ensimmäisen alueen
 * alialueeksi toiseen alueeseen.
 * \param id ensimmäinen alue
 * \param parentid toinen alue
 * \return totuusarvo true tai false
 */
bool Datastructures::add_subregion_to_region(RegionID id, RegionID parentid)
{
    Region* & A = regionsdata_.at(id);
    Region* & B = regionsdata_.at(parentid);

    if((regionsdata_.find(id) == regionsdata_.end() or
        regionsdata_.find(parentid) == regionsdata_.end()) or
            (regionsdata_[id]->parent != nullptr))
        return false;
    else
    {
        regionsdata_[id]->parent=B;
        regionsdata_[parentid]->children.insert(A);
        return true;
    }
}

/*!
 * \brief Datastructures::add_station_to_region Lisää aseman annettuun
 * alueeseen.
 * \param id aseman id
 * \param parentid alueen id
 * \return totuusarvo true tai false
 */
bool Datastructures::add_station_to_region(StationID id, RegionID parentid)
{
    if((stationdata_.find(id) == stationdata_.end()) or
        (regionsdata_.find(parentid) == regionsdata_.end()) or
            (stationdata_[id]->region != NO_REGION))
        return false;

    else
    {
        stationdata_[id]->region=parentid;
        regionsdata_[parentid]->stations.push_back(id);
        return true;
    }
}

/*!
 * \brief Datastructures::find_regions Etsii reitin annetusta solmusta juureen
 * \param current annettu solmu
 * \param vec Vektori johon reitti tallennetaan
 * \return Vektori johon reitti tallennetaan
 */
std::vector<RegionID> Datastructures::find_regions(Region* current, std::vector<RegionID> vec)
{
    while(current->parent != nullptr)
    {
        vec.push_back(current->regionid);
        current = current->parent;
    }
    vec.push_back(current->regionid);
    return vec;
}

/*!
 * \brief Datastructures::station_in_regions Palauttaa kaikki alueet, joihin
 * annettu asema kuuluu suoraan tai epäsuorasti. Paluuarvossa on ensin alue,
 * johon annettu asema kuuluu suoraan, sitten alue, johon tämä alue kuuluu jne
 * \param id aseman id numero
 * \return Vektori johon talletettu alueiden id numerot kuvauksen mukaisessa
 * järjestyksessä.
 */
std::vector<RegionID> Datastructures::station_in_regions(StationID id)
{
    std::vector<RegionID> vec;

    if(stationdata_.find(id) != stationdata_.end())
    {
        if(stationdata_[id]->region == NO_REGION)
            return {};
        else
        {

            Region* current = regionsdata_[stationdata_.at(id)->region];

            return find_regions(current,vec);
        }
    }
    else
        return {NO_REGION};
}

/*!
 * \brief Datastructures::list_subregions Etsii ja tallentaa reitin alueesta juureen.
 * \param region Alueen struct
 * \param reglist vektori johon alueiden id numerot tallennetaan kuvaillussa
 * järjestyksessä.
 */
void Datastructures::list_subregions(Region* region, std::vector<RegionID>& reglist)
{
    std::set<Region*> :: iterator reg;
    for(reg = region->children.begin(); reg != region->children.end(); reg++)
    {
        reglist.push_back((*reg)->regionid);
        list_subregions(*reg, reglist);
    }

}
/*!
 * \brief Datastructures::all_subregions_of_region Palauttaa kaikki alueet,
 * jotka kuuluvat annettuun alueeseen suoraan tai epäsuorasti alialueina
 * \param id alueen id numero
 * \return Vektori johon talletettu alueiden id numerot kuvauksen mukaisessa
 * järjestyksessä.
 */
std::vector<RegionID> Datastructures::all_subregions_of_region(RegionID id)
{
    std::vector<RegionID> subreg;


    if(regionsdata_.find(id) != regionsdata_.end())
    {
        if(regionsdata_[id]->children.empty())
            return {};

        else
        {
            list_subregions(regionsdata_[id], subreg);
            return subreg;
        }


    }
    else
        return {NO_REGION};
}

/*!
 * \brief Datastructures::stations_closest_to Palauttaa etäisyysjärjestyksessä
 * kolme annettua koordinaattia lähinnä olevaa asemaa.
 * \param xy annetut koordinaatit
 * \return Vektorin jossa ehdot täyttävät asemat
 */
std::vector<StationID> Datastructures::stations_closest_to(Coord xy)
{
    std::multimap<Distance, StationID> closeststations;
    std::vector<StationID> closestids;

    for(auto & entry: stationdata_)
    {
        int distx = (entry.second->coord.x - xy.x);
        int disty = (entry.second->y_coord - xy.y);
        closeststations.insert({calculate_distance(distx, disty), entry.first});
        if(closeststations.size() > 3)
            closeststations.erase(prev(closeststations.end()));
    }

    for(auto & entry: closeststations)
        closestids.push_back(entry.second);


    return closestids;
}

/*!
 * \brief Datastructures::remove_station Poistaa annetulla id:llä olevan
 * aseman, mikäli moinen löytyy
 * \param id aseman id numero
 * \return totuusarvo true tai false
 */
bool Datastructures::remove_station(StationID id)
{
    if(stationdata_.find(id) != stationdata_.end())
    {
        sortingstation_.erase(stationdata_[id]->name);

        auto station = stationdata_.find(id);
        delete (station->second);
        stationdata_.erase(id);

        return true;
    }

    else
        return false;
}

/*!
 * \brief Datastructures::common_parent_of_regions Palauttaa ”lähimmän” alueen
 * aluehierarkiassa, jonka alialueita molemmat annetut alueet ovat suoraan tai
 * epäsuorasti.
 * \param id1 annettu alue 1
 * \param id2annettu alue 2
 * \return Lähin alue joissa molemmat ovat alialueina, mikäli moista ei ole
 * plautetaan NO_COOR.
 */
RegionID Datastructures::common_parent_of_regions(RegionID id1, RegionID id2)
{
    std::vector<RegionID> path1;
    std::vector<RegionID> path2;
    std::vector<RegionID> commonpath;


    if(regionsdata_.find(id1) != regionsdata_.end() and regionsdata_.find(id2) != regionsdata_.end())
    {
        Region* current1 = regionsdata_[id1];

        path1 = find_regions(current1,path1);

        Region* current2 = regionsdata_[id2];

        path2 = find_regions(current2,path2);

        auto it1 = path1.size()-1;
        auto it2 = path2.size()-1;

        while(path1[it1] == path2[it2] and it1 >0 and it2 >0)
        {
            commonpath.push_back(path1[it1]);
            it1--;
            it2--;
        }
        if(commonpath.empty())
            return NO_REGION;

        else
            return commonpath.back();
    }
    else
        return NO_REGION;
}

/*!
 * \brief Datastructures::add_train Lisää tietorakenteeseen uuden junan annetulla uniikilla
 * id:llä. Juna kulkee annettujen asemien läpi ja lähtee niiltä annettuina
 * aikoina. Viimeisen aseman aika on saapumisaika pääteasemalle
 * \param trainid junan id
 * \param stationtimes aikataulut
 * \return true or false
 */
bool Datastructures::add_train(TrainID trainid, std::vector<std::pair<StationID, Time> > stationtimes)
{
    std::vector<std::pair<StationID, Station*>> stations;
    std::vector<StationID> list_of_stationids;

    if(traindata_.find(trainid) == traindata_.end())
    {
        for(auto & entry : stationtimes)
        {

            if(stationdata_.find(entry.first) != stationdata_.end())
            {
                stations.push_back({entry.first, stationdata_.at(entry.first)});
                list_of_stationids.push_back(entry.first);

            }
            else
                return false;
        }

        Train* new_train = new Train{trainid, stationtimes, stations, list_of_stationids};
        traindata_.insert({trainid, new_train});

        for(auto & entry2 : stationtimes)
        {
            stationdata_.at(entry2.first)->trains.insert({trainid, new_train});
            add_departure(entry2.first, trainid, entry2.second);
        }

        return true;
    }

    return false;
}

/*!
 * \brief Datastructures::next_stations_from Palauttaa asemat, jotka ovat heti
 * seuraavia asemia annetun aseman läpi kulkevissa junayhteyksissä.
 * \param id aseman id
 * \return  vektorin jossa seuraavien asemien idt
 */
std::vector<StationID> Datastructures::next_stations_from(StationID id)
{

    std::vector<StationID> next_stations;
    std::vector<StationID> routevector;

    if(stationdata_.find(id) != stationdata_.end())
    {
        for(auto & entry : stationdata_[id]->trains)
        {
            for(unsigned long i = 0; i < entry.second->stationids.size() -1; i++)
                if(entry.second->stationids[i] == id)
                    next_stations.push_back(entry.second->stationids[i+1]);

        }

        return next_stations;
    }

    return {NO_STATION};

}

/*!
 * \brief Datastructures::train_stations_from Palauttaa luettelon asemista,
 * joiden läpi annettu juna kulkee lähdettyään annetulta asemalta.
 * \param stationid aseman id
 * \param trainid junan id
 * \return vektori junan asemista
 */
std::vector<StationID> Datastructures::train_stations_from(StationID stationid, TrainID trainid)
{
    std::vector<StationID> stations_after;

    if(stationdata_.find(stationid) != stationdata_.end() and
            (stationdata_[stationid]->trains.find(trainid)
             != stationdata_[stationid]->trains.end()))
    {
        stations_after = traindata_[trainid]->stationids;

        int i = 0;
        while(stations_after[i] != stationid)
        stations_after.erase(stations_after.begin());

        stations_after.erase(stations_after.begin());
        if(!stations_after.empty())
            return stations_after;
        else
            return {NO_STATION};
}

    return {NO_STATION};
}

/*!
 * \brief Datastructures::clear_trains tyhjentää juna tietorakenteen
 */
void Datastructures::clear_trains()
{

    {
    for(std::unordered_map<std::string, Train*>::iterator train = traindata_.begin();
        train != traindata_.end(); train++)
        {
        delete (train->second);
        }
    traindata_.clear();
    }

    for(auto & entry: stationdata_)
        {
        entry.second->depatures={};
        entry.second->trains={};

        }
}

/*!
 * \brief Datastructures::route_any Palauttaa jonkin (mielivaltaisen) reitin
 * annettujen asemien välillä
 * \param fromid lähtöasema
 * \param toid määränpää
 * \return vektori jossa asemien idt ja välinen matka
 */
std::vector<std::pair<StationID, Distance>> Datastructures::route_any(StationID fromid, StationID toid)
{
    std::vector<std::pair<StationID, Distance>> any_route;
    std::stack<Station*> route_from_to;

    // Tarkastetaan että asmemat ovat olemassa
    if(stationdata_.find(fromid) != stationdata_.end() and
            stationdata_.find(toid) != stationdata_.end())
    {
        // Asetetaan kaikkien asemien väriksi white
        for(auto it = stationdata_.begin(); it != stationdata_.end(); it++)
        {
            it->second->colour="white";
            it->second->from_node=nullptr;
            it->second->distance_between=0;
            it->second->neighbours={};
        }


        // Asetetaan lähtöaseman etäisyys 0 Ja haetaan naapuriasemat
        // next_station_from funktion avulla
        stationdata_[fromid]->distance_between = 0;
        stationdata_[fromid]->neighbours = next_stations_from(fromid);

        // Lisätään lähtöasema pinoon sekä palautus vektoriin
        route_from_to.push(stationdata_[fromid]);
        //any_route.push_back({fromid, stationdata_[fromid]->distance_between});

        // Nimetään pinon päällimmäinen u:ksi
        auto u = route_from_to.top();

        // Kunnes pino on tyjä tai annettu määränpää saavutetaan
        while(!route_from_to.empty() and u!=stationdata_[toid])
        {

            // Päällimmäinen alkio mikäli white
            if(u->colour=="white")
            {
                // vaihdetaan väri ja lisätään u:n naapurit pinoon
                u->colour="grey";

                auto v = u->neighbours;

                for (auto & element : v)
                {
                    // Tallennetaan etäisyys, naapurit sekä tuloasema
                    // jokaiselle u:n naapurille
                    stationdata_[element]->distance_between = u->distance_between + calculate_distance(
                                    (stationdata_[element]->X_coord -  u->X_coord),
                                    (stationdata_[element]->y_coord -  u->y_coord));
                    stationdata_[element]->from_node=u;
                    stationdata_[element]->neighbours=next_stations_from(element);

                    // Lisätään u:n naapurit pinoon sekä palautusvektoriin
                    // mikäli naapuri ei siellä  vielä ole                  
                    route_from_to.push(stationdata_[element]);

                    }
                u=route_from_to.top();
                }
            // Mikäli u on jo kertaalleen käsitelty poistetaan u pinosta
            else if(u->colour=="grey")
            {
                route_from_to.pop();
            }
            // Mikäli väri ei ole white tai grey palautetaan tyhjä vektori
            else
            {
                return {};
            }
        }

        // Mikäli reittiä ei löydy palautetaan tyhjä vektori,
        // muuten palautusvektori
        if(u->stationid!=toid)
            return {};

        else
        {
            std::vector<std::pair<StationID, Distance>> any_route2;

            while(u!=stationdata_[fromid])
            {
                any_route.push_back({u->stationid, u->distance_between});
                u=u->from_node;
            }
            any_route.push_back({stationdata_[fromid]->stationid,
                                 stationdata_[fromid]->distance_between});

            for(auto i=any_route.rbegin(); i<any_route.rend(); i++)
                any_route2.push_back(*i);

            return any_route2;
        }
    }
    // Mikäli jompaakumpaa annettua asemaa ei löydy
    else
    {
        any_route.push_back({NO_STATION, NO_DISTANCE});
        return any_route;
    }
}

/*!
 * \brief Datastructures::route_least_stations Palauttaa sellaisen reitin
 * annettujen asemien välillä, jossa on mahdollisimman vähän asemia
 * \param fromid lähtöasema
 * \param toid määränpää
 * \return reitin asemien idt ja etäisyydet
 */
std::vector<std::pair<StationID, Distance>> Datastructures::route_least_stations(StationID fromid, StationID toid)
{
    std::vector<std::pair<StationID, Distance>> shortest_route;
    std::list<Station*> shortest_route_from_to;

    // Tarkastetaan että asmemat ovat olemassa
    if(stationdata_.find(fromid) != stationdata_.end() and
            stationdata_.find(toid) != stationdata_.end())
    {
        for(auto it = stationdata_.begin(); it != stationdata_.end(); it++)
        {
            it->second->colour="white";
            it->second->from_node=nullptr;
            it->second->distance_between=0;
            it->second->neighbours={};
        }

        stationdata_[fromid]->colour="grey";
        stationdata_[fromid]->neighbours = next_stations_from(fromid);
        shortest_route_from_to.push_back(stationdata_[fromid]);

        auto u = shortest_route_from_to.front();

        while (!shortest_route_from_to.empty() and u->stationid!=toid)
        {
            shortest_route_from_to.pop_front();
            u->colour="grey";
            auto v = u->neighbours;

            for (auto & element : v)
            {

                if(stationdata_[element]->stationid==toid)
                {
                    {
                        std::vector<std::pair<StationID, Distance>> any_route2;

                        stationdata_[element]->distance_between = u->distance_between + calculate_distance(
                                        (stationdata_[element]->X_coord -  u->X_coord),
                                        (stationdata_[element]->y_coord -  u->y_coord));
                        stationdata_[element]->from_node=u;

                        shortest_route_from_to.push_back(stationdata_[element]);

                        while(u!=stationdata_[fromid])
                        {
                            shortest_route.push_back({u->stationid, u->distance_between});
                            u=u->from_node;
                        }
                        shortest_route.push_back({stationdata_[fromid]->stationid,
                                             stationdata_[fromid]->distance_between});

                        for(auto i=shortest_route.rbegin(); i<shortest_route.rend(); i++)
                            any_route2.push_back(*i);

                        any_route2.push_back({toid, stationdata_[toid]->distance_between});

                        return any_route2;
                    }
                }
                else if(stationdata_[element]->colour=="white")
                {
                stationdata_[element]->distance_between = u->distance_between + calculate_distance(
                                (stationdata_[element]->X_coord -  u->X_coord),
                                (stationdata_[element]->y_coord -  u->y_coord));
                stationdata_[element]->from_node=u;
                stationdata_[element]->neighbours=next_stations_from(element);
                stationdata_[element]->colour="grey";

                    // Lisätään u:n naapurit pinoon sekä palautusvektoriin
                    // mikäli naapuri ei siellä  vielä ole
                shortest_route_from_to.push_back(stationdata_[element]);

                }
                else if(stationdata_[element]->colour=="grey")
                    stationdata_[element]->colour="black";


            }

                u=shortest_route_from_to.front();

            }

        return {};

    }
    // Mikäli jompaakumpaa annettua asemaa ei löydy
    else
    {
        shortest_route.push_back({NO_STATION, NO_DISTANCE});
        return shortest_route;
    }
}

/*!
 * \brief Datastructures::route_with_cycle Palauttaa annetulta asemalta
 * lähtevän reitin, jossa on sykli, ts. reitti päätyy uudelleen jollekin
 * reitillä jo oleva asemalle.
 * \param fromid lähtöasema
 * \return vektorin jossa syklin asemien idt
 */
std::vector<StationID> Datastructures::route_with_cycle(StationID fromid)
{
    std::vector<StationID> any_route;
    std::stack<Station*> route_from_to;


    // Tarkastetaan että asmema on olemassa
    if(stationdata_.find(fromid) != stationdata_.end() )
    {
        // Nollataan asemat
        for(auto it = stationdata_.begin(); it != stationdata_.end(); it++)
        {
            it->second->colour="white";
            it->second->from_node=nullptr;
            it->second->distance_between=0;
            it->second->neighbours={};
        }

        // Haetaan naapuriasemat next_station_from funktion avulla
        stationdata_[fromid]->neighbours = next_stations_from(fromid);

        // Lisätään lähtöasema pinoon
        route_from_to.push(stationdata_[fromid]);

        // Nimetään pinon päällimmäinen u:ksi
        auto u = route_from_to.top();

        // Kunnes pino on tyjä tai silmukka löytyy
        while(!route_from_to.empty())
        {

            // Päällimmäinen alkio mikäli white
            if(u->colour=="white")
            {
                // vaihdetaan väri ja lisätään u:n naapurit pinoon
                u->colour="grey";

                auto v = u->neighbours;

                for (auto & element : v)
                {
                    // Mikäli naapuri white
                    if(stationdata_[element]->colour=="white")
                    {
                        // Tallennetaan naapurit sekä tuloasema jokaiselle u:n naapurille
                        stationdata_[element]->from_node=u;
                        stationdata_[element]->neighbours=next_stations_from(element);

                        // Lisätään u:n naapurit pinoon
                        route_from_to.push(stationdata_[element]);
                    }
                    else
                    {
                        std::vector<StationID> any_route2;

                        while(u!=stationdata_[fromid])
                        {
                            any_route.push_back(u->stationid);
                            u=u->from_node;
                        }
                        any_route.push_back(stationdata_[fromid]->stationid);

                        for(auto i=any_route.rbegin(); i<any_route.rend(); i++)
                            any_route2.push_back(*i);

                        any_route2.push_back(stationdata_[element]->stationid);

                        return any_route2;

                    }

                    }

                u=route_from_to.top();
                }

            // Mikäli väri ei ole white tai grey palautetaan tyhjä vektori
            else
            {
                route_from_to.pop();
            }
        }

        // Mikäli reittiä ei löydy palautetaan tyhjä vektori,
        // muuten palautusvektori
        return {};

    }
    // Mikäli jompaakumpaa annettua asemaa ei löydy
    else
    {
        any_route.push_back({NO_STATION});
        return any_route;
    }

}


std::vector<std::pair<StationID, Distance>> Datastructures::route_shortest_distance(StationID fromid, StationID toid)
{
    std::vector<std::pair<StationID, Distance>> shortest_route;
    std::list<Station*> shortest_route_from_to;

    // Tarkastetaan että asmemat ovat olemassa
    if(stationdata_.find(fromid) != stationdata_.end() and
            stationdata_.find(toid) != stationdata_.end())
    {
        for(auto it = stationdata_.begin(); it != stationdata_.end(); it++)
        {
            it->second->colour="white";
            it->second->from_node=nullptr;
            it->second->distance_between=0;
            it->second->neighbours={};
        }

        stationdata_[fromid]->colour="grey";
        stationdata_[fromid]->neighbours = next_stations_from(fromid);
        shortest_route_from_to.push_back(stationdata_[fromid]);

        auto u = shortest_route_from_to.front();

        while (!shortest_route_from_to.empty())
        {
            shortest_route_from_to.pop_front();
            u->colour="grey";
            auto v = u->neighbours;

            for (auto & element : v)
            {

                if(stationdata_[element]->stationid==toid)
                {
                    {
                        std::vector<std::pair<StationID, Distance>> any_route2;

                        stationdata_[element]->distance_between = u->distance_between + calculate_distance(
                                        (stationdata_[element]->X_coord -  u->X_coord),
                                        (stationdata_[element]->y_coord -  u->y_coord));
                        stationdata_[element]->from_node=u;

                        shortest_route_from_to.push_back(stationdata_[element]);

                        while(u!=stationdata_[fromid])
                        {
                            shortest_route.push_back({u->stationid, u->distance_between});
                            u=u->from_node;
                        }
                        shortest_route.push_back({stationdata_[fromid]->stationid,
                                             stationdata_[fromid]->distance_between});

                        for(auto i=shortest_route.rbegin(); i<shortest_route.rend(); i++)
                            any_route2.push_back(*i);

                        any_route2.push_back({toid, stationdata_[toid]->distance_between});

                        return any_route2;
                    }
                }
                else if(stationdata_[element]->colour=="white")
                {
                stationdata_[element]->distance_between = u->distance_between + calculate_distance(
                                (stationdata_[element]->X_coord -  u->X_coord),
                                (stationdata_[element]->y_coord -  u->y_coord));
                stationdata_[element]->from_node=u;
                stationdata_[element]->neighbours=next_stations_from(element);
                stationdata_[element]->colour="grey";

                    // Lisätään u:n naapurit pinoon sekä palautusvektoriin
                    // mikäli naapuri ei siellä  vielä ole
                shortest_route_from_to.push_back(stationdata_[element]);

                }
                else if(stationdata_[element]->colour=="grey")
                    stationdata_[element]->colour="black";


            }

                u=shortest_route_from_to.front();

            }

        return {};

    }
    // Mikäli jompaakumpaa annettua asemaa ei löydy
    else
    {
        shortest_route.push_back({NO_STATION, NO_DISTANCE});
        return shortest_route;
    }
}
std::vector<std::pair<StationID, Time>> Datastructures::route_earliest_arrival(StationID /*fromid*/, StationID /*toid*/, Time /*starttime*/)
{
    // Replace the line below with your implementation
    // Also uncomment parameters ( /* param */ -> param )
    throw NotImplemented("route_earliest_arrival()");
}
