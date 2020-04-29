function fetchAPI(apiEndPoint) {
  fetch(scriptRoot + apiEndPoint)
    .then(response => {
      return response.json();
    })
    .then(data => {
      document.getElementById("json").textContent = JSON.stringify(data, undefined, 2);
      document.getElementById("results").classList.remove("invisible");
      document.getElementById("loader").classList.add("invisible");
    });
  document.getElementById("results").classList.add("invisible");
  document.getElementById("loader").classList.remove("invisible");

}

let search = () => {
  let m =  document.getElementById('search-input').value;
  fetchAPI(scriptRoot + "/api/getDetailsForSeq?seq=" + encodeURIComponent(m));
}

let fetchSEQBySpecimen = () => {
  fetchAPI("/api/getSEQCountbySpecimenSource");
}

let fetchSEQByLocation = () => {
  fetchAPI("/api/getSEQCountbyLocation");
}

let fetchSEQByTech = () => {
  fetchAPI("/api/getSEQCountbytech");
}

let fetchAllaccessions = () => {
  fetchAPI("/api/getAllaccessions");
};

/**
 * Show form if checked
 */
let fillFormSpot = document.getElementById('metadata_fill_form_spot');
function displayForm() {
  if (document.getElementById('metadata_form').checked) {
    fillFormSpot.classList.remove("invisible");
    return;
  }
  fillFormSpot.classList.add("invisible");
}
