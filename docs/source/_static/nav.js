/*
 * Mobile navigation bar.
 *
 * Alabaster has no mobile navigation. Below its 1180px breakpoint it drops the
 * whole sidebar in as a block above the content, which on a phone is roughly
 * 1500px of logo, description, 35 nav rows and a search box before the first
 * word of the page.
 *
 * This adds a sticky bar with a menu button and collapses the sidebar behind
 * it. The bar and the collapsing are both scoped to the same breakpoint in
 * style.css, so nothing here affects the desktop layout. If this script does
 * not run the stylesheet caps the sidebar height instead, so the page stays
 * usable without javascript.
 */
(function () {
  'use strict';

  var BREAKPOINT = 1180;

  function build() {
    var document_ = document.querySelector('div.document');
    var sidebar = document.querySelector('div.sphinxsidebar');
    if (!document_ || !sidebar || document.querySelector('.ms-navbar')) {
      return;
    }

    var title = (document.querySelector('div.sphinxsidebar h1 a') || {}).textContent
      || document.title.split('—').pop().trim()
      || 'mofstructure';

    var bar = document.createElement('nav');
    bar.className = 'ms-navbar';
    bar.setAttribute('aria-label', 'Site navigation');

    var brand = document.createElement('a');
    brand.className = 'ms-navbar-brand';
    brand.href = document.querySelector('link[rel="index"]') ? 'index.html' : '#';
    brand.textContent = title.trim();

    var button = document.createElement('button');
    button.className = 'ms-navbar-toggle';
    button.type = 'button';
    button.setAttribute('aria-expanded', 'false');
    button.setAttribute('aria-controls', 'ms-sidebar');
    button.setAttribute('aria-label', 'Open navigation menu');
    button.innerHTML = '<span class="ms-navbar-bars" aria-hidden="true"></span>';

    bar.appendChild(brand);
    bar.appendChild(button);
    document_.parentNode.insertBefore(bar, document_);

    sidebar.id = 'ms-sidebar';
    sidebar.classList.add('ms-collapsible');

    function setOpen(open) {
      sidebar.classList.toggle('ms-open', open);
      button.setAttribute('aria-expanded', open ? 'true' : 'false');
      button.setAttribute('aria-label',
        open ? 'Close navigation menu' : 'Open navigation menu');
      document.body.classList.toggle('ms-nav-open', open);
    }

    button.addEventListener('click', function () {
      setOpen(!sidebar.classList.contains('ms-open'));
    });

    // follow a link, close the menu behind you
    sidebar.addEventListener('click', function (event) {
      if (event.target.closest('a')) {
        setOpen(false);
      }
    });

    document.addEventListener('keydown', function (event) {
      if (event.key === 'Escape' && sidebar.classList.contains('ms-open')) {
        setOpen(false);
        button.focus();
      }
    });

    // widening past the breakpoint returns to the desktop layout
    window.addEventListener('resize', function () {
      if (window.innerWidth > BREAKPOINT) {
        setOpen(false);
      }
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', build);
  } else {
    build();
  }
})();
