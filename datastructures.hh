// Datastructures.hh
//
// Student name: Mikko Kaukonen


#ifndef DATASTRUCTURES_HH
#define DATASTRUCTURES_HH

#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <limits>
#include <functional>
#include <exception>
#include <map>
#include <set>
#include <unordered_set>

// Types for IDs
using StationID = std::string;
using TrainID = std::string;
using RegionID = unsigned long long int;
using Name = std::string;
using Time = unsigned short int;

// Return values for cases where required thing was not found
StationID const NO_STATION = "---";
TrainID const NO_TRAIN = "---";
RegionID const NO_REGION = -1;
Name const NO_NAME = "!NO_NAME!";
Time const NO_TIME = 9999;

// Return value for cases where integer values were not found
int const NO_VALUE = std::numeric_limits<int>::min();


// Type for a coordinate (x, y)
struct Coord
{
    int x = NO_VALUE;
    int y = NO_VALUE;
};

// Example: Defining == and hash function for Coord so that it can be used
// as key for std::unordered_map/set, if needed
inline bool operator==(Coord c1, Coord c2) { return c1.x == c2.x && c1.y == c2.y; }
inline bool operator!=(Coord c1, Coord c2) { return !(c1==c2); } // Not strictly necessary

struct CoordHash
{
    std::size_t operator()(Coord xy) const
    {
        auto hasher = std::hash<int>();
        auto xhash = hasher(xy.x);
        auto yhash = hasher(xy.y);
        // Combine hash values (magic!)
        return xhash ^ (yhash + 0x9e3779b9 + (xhash << 6) + (xhash >> 2));
    }
};

// Example: Defining < for Coord so that it can be used
// as key for std::map/set
inline bool operator<(Coord c1, Coord c2)
{
    if (c1.y < c2.y) { return true; }
    else if (c2.y < c1.y) { return false; }
    else { return c1.x < c2.x; }
}

// Return value for cases where coordinates were not found
Coord const NO_COORD = {NO_VALUE, NO_VALUE};

// Type for a distance (in metres)
using Distance = int;

// Return value for cases where Distance is unknown
Distance const NO_DISTANCE = NO_VALUE;

// This exception class is there just so that the user interface can notify
// about operations which are not (yet) implemented
class NotImplemented : public std::exception
{
public:
    NotImplemented() : msg_{} {}
    explicit NotImplemented(std::string const& msg) : msg_{msg + " not implemented"} {}

    virtual const char* what() const noexcept override
    {
        return msg_.c_str();
    }
private:
    std::string msg_;
};


// This is the class you are supposed to implement

class Datastructures
{
public:
    Datastructures();
    ~Datastructures();

    //
    // Old assignment 1 operations
    //


    unsigned int station_count();

    void clear_all();

    std::vector<StationID> all_stations();

    bool add_station(StationID id, Name const& name, Coord xy);

    Name get_station_name(StationID id);

    Coord get_station_coordinates(StationID id);

    // We recommend you implement the operations below only after implementing the ones above

    std::vector<StationID> stations_alphabetically();

    std::vector<StationID> stations_distance_increasing();

    StationID find_station_with_coord(Coord xy);

    bool change_station_coord(StationID id, Coord newcoord);

    bool add_departure(StationID stationid, TrainID trainid, Time time);

    bool remove_departure(StationID stationid, TrainID trainid, Time time);

    std::vector<std::pair<Time, TrainID>> station_departures_after(StationID stationid, Time time);

    // We recommend you implement the operations below only after implementing the ones above

    bool add_region(RegionID id, Name const& name, std::vector<Coord> coords);

    std::vector<RegionID> all_regions();

    Name get_region_name(RegionID id);

    std::vector<Coord> get_region_coords(RegionID id);

    bool add_subregion_to_region(RegionID id, RegionID parentid);

    bool add_station_to_region(StationID id, RegionID parentid);

    std::vector<RegionID> station_in_regions(StationID id);

    // Non-compulsory operations

    std::vector<RegionID> all_subregions_of_region(RegionID id);

    std::vector<StationID> stations_closest_to(Coord xy);

    bool remove_station(StationID id);

    RegionID common_parent_of_regions(RegionID id1, RegionID id2);

    //
    // New assignment 2 operations
    //

    bool add_train(TrainID trainid, std::vector<std::pair<StationID, Time>> stationtimes);

    std::vector<StationID> next_stations_from(StationID id);

    std::vector<StationID> train_stations_from(StationID stationid, TrainID trainid);

    void clear_trains();

    std::vector<std::pair<StationID, Distance>> route_any(StationID fromid, StationID toid);

    // Non-compulsory operations

    std::vector<std::pair<StationID, Distance>> route_least_stations(StationID fromid, StationID toid);

    std::vector<StationID> route_with_cycle(StationID fromid);

    std::vector<std::pair<StationID, Distance>> route_shortest_distance(StationID fromid, StationID toid);
    std::vector<std::pair<StationID, Time>> route_earliest_arrival(StationID fromid, StationID toid, Time starttime);

private:

    // Add stuff needed for your class implementation here
    struct Train;

    struct Station{

        StationID stationid;
        Name name;
        Coord coord;
        int Distance;
        int y_coord;
        int X_coord;
        std::map<Time,std::vector<TrainID>> depatures;
        RegionID region;
        std::map<TrainID, Train*> trains;

        int distance_between;
        Station* from_node;
        std::string colour;
        std::vector<StationID> neighbours;
        int distance_from_start;

    };

    struct Region{

        RegionID regionid;
        Name name;
        std::vector<Coord> coordinates;
        std::vector<StationID> stations;

        std::set<Region*> children;
        Region* parent; // nullptr if root
    };

    struct Train{

        TrainID trainid;
        std::vector<std::pair<StationID, Time>> route;
        std::vector<std::pair<StationID, Station*>> stations;
        std::vector<StationID> stationids;

    };

    std::unordered_map<StationID, Station*> stationdata_;
    std::multimap<Name, Station*> sortingstation_;
    std::unordered_map<unsigned long long int, Region*> regionsdata_;
    std::unordered_map<TrainID, Train*> traindata_;



    int calculate_distance(int x, int y);

    std::vector<RegionID> find_regions(Region* current, std::vector<RegionID> vec);

    void list_subregions(Region* regions, std::vector<RegionID>& reglist);

    std::vector<RegionID> path_to_root(RegionID);

};

#endif // DATASTRUCTURES_HH
