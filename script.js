$(document).ready(function() {
  // Configure MathJax.
  MathJax.Hub.Config({
    extensions: ["tex2jax.js"],
    jax: ["input/TeX", "output/HTML-CSS"],
    tex2jax: {
      inlineMath: [ ['$$','$$'], ["\\(","\\)"] ],
      displayMath: [ ["\\[","\\]"] ],
    },
    "HTML-CSS": { availableFonts: ["TeX"] }
  });

  // Configure Google+ button.
  var po = document.createElement('script');
  po.type = 'text/javascript';
  po.async = true;
  po.src = 'https://apis.google.com/js/plusone.js';
  var s = document.getElementsByTagName('script')[0];
  s.parentNode.insertBefore(po, s);
});
