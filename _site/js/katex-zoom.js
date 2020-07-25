function on_click_katex() {
  var elem = this.parentElement;
  
  var pack = document.createElement('span');
  var abso = document.createElement('span');
  var over = document.createElement('span');
  var zoom = document.createElement('span');
  var tag, tag2;
  
  $(pack).addClass('katex-zoom-frame');
  $(abso).addClass('katex-zoom-absolute');
  $(zoom).addClass('katex-zoom-main');
  $(over).addClass('katex-zoom-overlayer');
  $(elem).addClass('katex-zoom-tmp');

  $(abso).text('.').append(zoom);
  $(zoom).append(this.cloneNode(true));
  $(pack).append(over, abso);
  $(elem).prepend(pack);

  $(zoom).css('left', (this.offsetWidth - zoom.offsetWidth) / 2);
  $(zoom).css('top', (this.offsetHeight - zoom.offsetHeight) / 2);
  $(zoom).attr('id', 'katex-zoom');

  if ($('#katex-zoom .tag').length > 0) {
    $('#katex-zoom .tag')[0].remove();
    var tag = $('.katex-zoom-tmp .tag')[0];
    tag2 = tag.cloneNode(true);
    $(tag) .addClass('katex-zoom-tag');
    $(tag2).addClass('katex-zoom-tag2 katex');
    $(elem.parentElement).prepend(tag2);
  }

  $(over).on('click', function() {
    $(pack).remove();
    $(elem).removeClass('katex-zoom-tmp');
    if (typeof tag !== "undefined") {
      $(tag) .removeClass('katex-zoom-tag');
      $(tag2).remove();
    }
  });
}

function enlarge_katex_loaded() {
  $('span.katex').on('dblclick', on_click_katex);
}
