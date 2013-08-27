$(document).ready(function() {
  // Initialize top bar.
  var threshold = $('#top').position().top + $('#top').outerHeight();
  $('#topbar').css('display', 'table');
  $(window).scroll(function() {
    if ($(window).scrollTop() > threshold) {
      $('#topbar').slideDown('fast');
    } else {
      $('#topbar').slideUp('fast');
    }
  });
  $('#topbar').hide();
});
