
var map = L.map( 'mapid', {
    center: [51.505, -0.09],  // Default to U.S.A
    minZoom: 2,
    zoom: 0
});

L.tileLayer( 'http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> | <a href="http://covid19.genenetwork.org/">COVID-19 PubSeq</a>',
    subdomains: ['a','b','c']
}).addTo(map);


function drawMap(){
    var mymap = map;

    fetch(scriptRoot + "api/getCountByGPS")
        .then(response => {
            console.log(response)
            return response.json();
        })
        .then(data => {
            /*
            for (var i=0; i<data.length;i++) {
                gps=data[i]["GPS"].split(" ")
                var circle = L.circle([gps[1], gps[0]], {
                    color: 'red',
                    fillColor: '#f03',
                    fillOpacity: 0.5,
                    radius: parseInt(data[i]["count"])  //not working for whatever reason
                }).addTo(mymap);
            }
            */
            updateMapMarkers(data);

      });
    document.getElementById("map_view").classList.remove("invisible");
    map.invalidateSize();
}



/* This function updates the map with markers
 *
*/
function updateMapMarkers(data) {
    let markers = L.markerClusterGroup();
    for (let i = 0; i < data.length; i++) {
        let {"count": fastaCount, GPS, LocationLabel: label } = data[i];
        let coordinates = GPS.split(" ");
        if (!(coordinates == null)) {
            let lat, lon;
            [lon, lat] = coordinates.map(parseFloat);
            let point = L.point()
            marker = (L.marker([lat, lon]));
            //  .bindPopup("<b>" + label + "</b><br/>" + "FastaCount: " +fastaCount));
            markers.addLayer(marker);
        }
    }
    map.addLayer(markers);
}
