fetch(scriptRoot + "/api/getAllaccessions")
  .then(response => {
    return response.json();
  })
  .then(data => {
    console.log('test');
    console.log(data);
  })

/**
 * Show form if checked
 */
let fillFormSpot = document.getElementById('metadata_fill_form_spot');
function displayForm() {
  if (document.getElementById('metadata_form').checked) {
    fillFormSpot.classList.remove("invisible");
  } else {
    fillFormSpot.classList.add("invisible");
    console.log("visible");
  }
}
