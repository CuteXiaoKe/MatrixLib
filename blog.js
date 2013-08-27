var posts = [
  ["Numerical analysis, linear algebra, and Java", "07/11/2011"],
  ["Numerical stability", "07/13/2011"],
  ["Gauss-Jordan reduction and determinants", "07/15/2011"],
  ["Stable matrix inversion", "07/17/2011"],
  ["Householder reflections", "07/19/2011"],
  ["Computing eigenvalues", "07/20/2011"],
  ["Computing eigenvectors", "07/24/2011"],
  ["Future directions for MatrixLib", "07/29/2011"],
  ["Precision in complex computation", "02/26/2012"],
];

$(document).ready(function() {
  $('#topbar').css('display', 'table');
  var id = window.location.hash.substr(1);
  if (id != "") {
    load_post(id);
  } else {
    load_post(posts.length);
  }
});

function load_post(id) {
  $('#blog-content').innerHTML = '';
  var name = "blog/" + id.toString() + ".html";
  $("#blog-content").load(name, function() {
    var title = '<p class="section">' + posts[id-1][0] + '</p>';
    var date = '<p class="date">' + posts[id-1][1] + '</p>';
    $("#blog-content").prepend(title + date);
    MathJax.Hub.Queue(["Typeset", MathJax.Hub, "blog-content"]);
  });
}

var showing = false;

function toggle_posts() {
  if (showing == true) {
    $('#posts').slideUp('fast');
    $('#posts').empty();
    showing = false;
  } else {
    for (var i = posts.length - 1; i >= 0; i--) {
      var title = '<span class="section">' + 
        '<a href="blog.html#' + (i+1) + '">' + posts[i][0] + '</a></span>';
      var date = '<span class="date">' + posts[i][1] + '</span>';
      $('#posts').append('<div>' + title + ' ' + date + '</div>');
      $('#posts').hide().slideDown('fast');
    }
    showing = true;
  }
}
